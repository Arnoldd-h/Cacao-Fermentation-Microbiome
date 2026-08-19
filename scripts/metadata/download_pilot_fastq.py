#!/usr/bin/env python3
"""Download the pilot FASTQ manifest resumably and validate size plus ENA MD5.

Inputs: metadata/pilot_manifest.tsv and config/config.yaml. Raw FASTQ are written
under data/raw/ and remain excluded from Git. Output:
results/qc/pilot_download_validation.tsv. Fails after recording any download,
disk-space, size or checksum error.
"""

from __future__ import annotations

import argparse
import logging
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_json_yaml  # noqa: E402
from cacao_inventory.download import (  # noqa: E402
    download_file,
    expand_pilot_manifest,
    file_md5,
)
from cacao_inventory.io import read_tsv, write_tsv_atomic  # noqa: E402
from cacao_inventory.schema import (  # noqa: E402
    DOWNLOAD_VALIDATION_COLUMNS,
    PILOT_MANIFEST_COLUMNS,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "config.yaml")
    parser.add_argument(
        "--manifest", type=Path, default=ROOT / "metadata" / "pilot_manifest.tsv"
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data" / "raw")
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "results" / "qc" / "pilot_download_validation.tsv",
    )
    return parser.parse_args()


def _git_commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
    )
    return completed.stdout.strip()


def _display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    config = load_json_yaml(args.config)
    settings = config.get("pilot_download")
    if not isinstance(settings, dict):
        raise ValueError("config.yaml must define pilot_download")
    chunk_size = int(settings["chunk_size_bytes"])
    free_space_factor = float(settings["minimum_free_space_factor"])
    if chunk_size <= 0 or free_space_factor < 1:
        raise ValueError("Invalid pilot_download settings")

    manifest = read_tsv(args.manifest, PILOT_MANIFEST_COLUMNS)
    records = expand_pilot_manifest(manifest)
    total_bytes = sum(int(record["expected_bytes"]) for record in records)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    free_bytes = shutil.disk_usage(args.output_dir).free
    required_bytes = int(total_bytes * free_space_factor)
    if free_bytes < required_bytes:
        raise OSError(
            f"Insufficient disk space: {free_bytes} free bytes, {required_bytes} required"
        )
    logging.info(
        "Pilot requires %d FASTQ files and %d bytes; %d bytes are free",
        len(records),
        total_bytes,
        free_bytes,
    )

    validated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    git_commit = _git_commit()
    python_version = platform.python_version()
    report_rows: list[dict[str, str]] = []
    failures: list[str] = []
    for record in records:
        destination = args.output_dir / record["target_relative"]
        status = "failed"
        observed_bytes = ""
        observed_md5 = ""
        try:
            logging.info("Downloading or validating %s", record["read_file"])
            status = download_file(
                record["source_url"],
                destination,
                int(record["expected_bytes"]),
                record["expected_md5"],
                chunk_size,
            )
            observed_bytes = str(destination.stat().st_size)
            observed_md5 = file_md5(destination, chunk_size)
        except Exception as exc:  # keep a complete resumable execution report
            failures.append(f"{record['read_file']}: {exc}")
            logging.error("%s", failures[-1])
        report_rows.append(
            {
                "study_id": record["study_id"],
                "bioproject": record["bioproject"],
                "sample_id": record["sample_id"],
                "run_accession": record["run_accession"],
                "read_file": record["read_file"],
                "source_url": record["source_url"],
                "local_path": _display_path(destination),
                "expected_bytes": record["expected_bytes"],
                "observed_bytes": observed_bytes,
                "expected_md5": record["expected_md5"],
                "observed_md5": observed_md5,
                "status": status,
                "validated_at_utc": validated_at,
                "git_commit": git_commit,
                "python_version": python_version,
            }
        )
    write_tsv_atomic(args.report, report_rows, DOWNLOAD_VALIDATION_COLUMNS)
    if failures:
        raise RuntimeError(f"{len(failures)} FASTQ files failed validation; see {args.report}")
    logging.info("Validated %d FASTQ files", len(report_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
