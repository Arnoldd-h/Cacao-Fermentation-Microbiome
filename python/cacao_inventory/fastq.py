"""Dependency-free gzip and four-line FASTQ integrity validation."""

from __future__ import annotations

import gzip
from pathlib import Path


def validate_fastq_file(path: str | Path) -> dict[str, int]:
    """Read a gzipped FASTQ fully and return structural summary metrics."""

    source = Path(path)
    read_count = 0
    base_count = 0
    minimum_length: int | None = None
    maximum_length = 0
    with gzip.open(source, "rb") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence_line = handle.readline()
            plus = handle.readline()
            quality_line = handle.readline()
            record_number = read_count + 1
            if not sequence_line or not plus or not quality_line:
                raise ValueError(f"Truncated FASTQ record {record_number} in {source}")
            sequence = sequence_line.rstrip(b"\r\n")
            quality = quality_line.rstrip(b"\r\n")
            if not header.startswith(b"@"):
                raise ValueError(f"Invalid FASTQ header at record {record_number} in {source}")
            if not plus.startswith(b"+"):
                raise ValueError(f"Invalid FASTQ separator at record {record_number} in {source}")
            if not sequence:
                raise ValueError(f"Empty sequence at record {record_number} in {source}")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"Sequence/quality length mismatch at record {record_number} in {source}"
                )
            length = len(sequence)
            read_count += 1
            base_count += length
            minimum_length = length if minimum_length is None else min(minimum_length, length)
            maximum_length = max(maximum_length, length)
    if read_count == 0:
        raise ValueError(f"FASTQ file contains no records: {source}")
    return {
        "read_count": read_count,
        "base_count": base_count,
        "minimum_read_length": minimum_length or 0,
        "maximum_read_length": maximum_length,
    }
