#!/usr/bin/env python3
"""Validate gzip integrity and four-line FASTQ structure for the pilot files.

Input: results/qc/pilot_download_validation.tsv. Output:
results/qc/pilot_fastq_validation.tsv. Every input must already have a valid
download status and resolve inside data/raw/. The validator reads every record,
thereby checking the gzip stream CRC, headers, separators and quality lengths.
"""

from __future__ import annotations

import argparse
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.fastq import validate_fastq_file  # noqa: E402
from cacao_inventory.io import read_tsv, write_tsv_atomic  # noqa: E402
from cacao_inventory.schema import (  # noqa: E402
    DOWNLOAD_VALIDATION_COLUMNS,
    FASTQ_VALIDATION_COLUMNS,
)


VALID_DOWNLOAD_STATUSES = {"downloaded_valid", "resumed_valid", "reused_valid"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--download-report",
        type=Path,
        default=ROOT / "results" / "qc" / "pilot_download_validation.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "qc" / "pilot_fastq_validation.tsv",
    )
    return parser.parse_args()


def _git_commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
    )
    return completed.stdout.strip()


def main() -> int:
    args = parse_args()
    download_rows = read_tsv(args.download_report, DOWNLOAD_VALIDATION_COLUMNS)
    raw_root = (ROOT / "data" / "raw").resolve()
    validated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    git_commit = _git_commit()
    python_version = platform.python_version()
    output_rows: list[dict[str, str]] = []
    failures: list[str] = []
    for row in download_rows:
        local_path = (ROOT / row["local_path"]).resolve()
        metrics = {
            "read_count": 0,
            "base_count": 0,
            "minimum_read_length": 0,
            "maximum_read_length": 0,
        }
        gzip_valid = "false"
        fastq_valid = "false"
        status = "failed"
        try:
            if row["status"] not in VALID_DOWNLOAD_STATUSES:
                raise ValueError(f"Download status is not valid: {row['status']}")
            if raw_root != local_path and raw_root not in local_path.parents:
                raise ValueError(f"FASTQ path is outside data/raw: {local_path}")
            metrics = validate_fastq_file(local_path)
            gzip_valid = "true"
            fastq_valid = "true"
            status = "valid"
        except Exception as exc:
            failures.append(f"{row['read_file']}: {exc}")
        output_rows.append(
            {
                "study_id": row["study_id"],
                "bioproject": row["bioproject"],
                "sample_id": row["sample_id"],
                "run_accession": row["run_accession"],
                "read_file": row["read_file"],
                "local_path": row["local_path"],
                "compressed_bytes": row["observed_bytes"],
                "read_count": str(metrics["read_count"]),
                "base_count": str(metrics["base_count"]),
                "minimum_read_length": str(metrics["minimum_read_length"]),
                "maximum_read_length": str(metrics["maximum_read_length"]),
                "gzip_valid": gzip_valid,
                "fastq_valid": fastq_valid,
                "status": status,
                "validated_at_utc": validated_at,
                "git_commit": git_commit,
                "python_version": python_version,
            }
        )
    write_tsv_atomic(args.output, output_rows, FASTQ_VALIDATION_COLUMNS)
    if failures:
        raise RuntimeError("; ".join(failures))
    print(f"Validated {len(output_rows)} gzip FASTQ files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
