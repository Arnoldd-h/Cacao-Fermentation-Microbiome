"""Resumable FASTQ download and checksum validation for a pilot manifest."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO, Callable
from urllib.parse import urlparse
from urllib.request import Request, urlopen


def _split(value: str) -> list[str]:
    return [part.strip() for part in value.split(";") if part.strip()]


def expand_pilot_manifest(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Expand paired manifest rows into one validated record per FASTQ file."""

    records: list[dict[str, str]] = []
    targets: set[str] = set()
    for row in rows:
        urls = _split(row["fastq_ftp"])
        checksums = _split(row["fastq_md5"])
        byte_values = _split(row["fastq_bytes"])
        if not urls or not (len(urls) == len(checksums) == len(byte_values)):
            raise ValueError(f"Unaligned FASTQ metadata for {row['run_accession']}")
        for url, checksum, byte_value in zip(urls, checksums, byte_values, strict=True):
            filename = PurePosixPath(urlparse(url).path).name
            if not filename:
                raise ValueError(f"FASTQ URL has no filename: {url}")
            try:
                expected_bytes = int(byte_value)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid FASTQ byte count for {row['run_accession']}: {byte_value}"
                ) from exc
            if expected_bytes <= 0:
                raise ValueError(f"FASTQ byte count must be positive for {row['run_accession']}")
            target = f"{row['study_id']}/{row['run_accession']}/{filename}"
            if target in targets:
                raise ValueError(f"Duplicate FASTQ target in manifest: {target}")
            targets.add(target)
            records.append(
                {
                    "study_id": row["study_id"],
                    "bioproject": row["bioproject"],
                    "sample_id": row["sample_id"],
                    "run_accession": row["run_accession"],
                    "read_file": filename,
                    "source_url": url,
                    "target_relative": target,
                    "expected_bytes": str(expected_bytes),
                    "expected_md5": checksum.lower(),
                }
            )
    return records


def file_md5(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    """Return the hexadecimal MD5 required by ENA FASTQ manifests."""

    digest = hashlib.md5(usedforsecurity=False)
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def validate_file(path: str | Path, expected_bytes: int, expected_md5: str) -> tuple[int, str]:
    """Fail clearly when a local file differs from its declared bytes or MD5."""

    source = Path(path)
    observed_bytes = source.stat().st_size
    if observed_bytes != expected_bytes:
        raise ValueError(
            f"Size mismatch for {source}: expected {expected_bytes}, observed {observed_bytes}"
        )
    observed_md5 = file_md5(source)
    if observed_md5.lower() != expected_md5.lower():
        raise ValueError(
            f"MD5 mismatch for {source}: expected {expected_md5}, observed {observed_md5}"
        )
    return observed_bytes, observed_md5


def download_file(
    url: str,
    destination: str | Path,
    expected_bytes: int,
    expected_md5: str,
    chunk_size: int,
    opener: Callable[[Request], BinaryIO] = urlopen,
) -> str:
    """Download through a .part file, resume with HTTP Range and validate before promotion."""

    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        validate_file(target, expected_bytes, expected_md5)
        return "reused_valid"

    partial = target.with_name(f"{target.name}.part")
    offset = partial.stat().st_size if partial.exists() else 0
    if offset > expected_bytes:
        raise ValueError(f"Partial file is larger than expected: {partial}")
    if offset == expected_bytes:
        validate_file(partial, expected_bytes, expected_md5)
        os.replace(partial, target)
        return "resumed_valid"

    request = Request(url, headers={"Range": f"bytes={offset}-"} if offset else {})
    with opener(request) as response:
        status = getattr(response, "status", None)
        append = bool(offset and status == 206)
        mode = "ab" if append else "wb"
        with partial.open(mode) as handle:
            while chunk := response.read(chunk_size):
                handle.write(chunk)

    validate_file(partial, expected_bytes, expected_md5)
    os.replace(partial, target)
    return "resumed_valid" if offset and append else "downloaded_valid"
