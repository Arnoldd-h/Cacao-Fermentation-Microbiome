#!/usr/bin/env python3
"""Build a validated machine-readable table from MultiQC FastQC output.

Inputs: MultiQC's multiqc_fastqc.txt and the pilot FASTQ validation report.
Output: results/qc/pilot/raw_read_quality.tsv by default. Fails if samples,
filenames or read counts disagree between FastQC and the validated FASTQ files.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.io import read_tsv, write_tsv_atomic  # noqa: E402
from cacao_inventory.provenance import qc_provenance
from cacao_inventory.qc import (  # noqa: E402
    MULTIQC_FASTQC_COLUMNS,
    build_raw_quality_rows,
)
from cacao_inventory.schema import (  # noqa: E402
    FASTQ_VALIDATION_COLUMNS,
    RAW_READ_QUALITY_COLUMNS,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--multiqc-fastqc",
        type=Path,
        default=ROOT
        / "results"
        / "qc"
        / "pilot"
        / "multiqc_raw"
        / "multiqc_report_data"
        / "multiqc_fastqc.txt",
    )
    parser.add_argument(
        "--fastq-validation",
        type=Path,
        default=ROOT / "results" / "qc" / "pilot_fastq_validation.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "qc" / "pilot" / "raw_read_quality.tsv",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    multiqc_rows = read_tsv(args.multiqc_fastqc, MULTIQC_FASTQC_COLUMNS)
    validation_rows = read_tsv(args.fastq_validation, FASTQ_VALIDATION_COLUMNS)
    output_rows = build_raw_quality_rows(multiqc_rows, validation_rows)
    write_tsv_atomic(args.output, output_rows, RAW_READ_QUALITY_COLUMNS)
    qc_provenance(
        args.output, root=ROOT,
        inputs=[args.multiqc_fastqc, args.fastq_validation, *(ROOT / row["local_path"] for row in validation_rows)],
        tools=("fastqc", "multiqc"), command=[sys.executable, *sys.argv],
    )
    logging.info("Recorded raw FastQC metrics for %d FASTQ files", len(output_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
