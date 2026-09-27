#!/usr/bin/env python3
"""Summarize post-Cutadapt FastQC and compare it with raw pilot QC.

Inputs: MultiQC FastQC data for interim FASTQ, the validated Cutadapt summary,
and raw_read_quality.tsv. Outputs: trimmed_read_quality.tsv and
read_quality_comparison.tsv. Fails when samples, directions, counts, paths, or
metadata disagree.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.io import read_tsv, write_tsv_atomic  # noqa: E402
from cacao_inventory.qc import (  # noqa: E402
    MULTIQC_FASTQC_COLUMNS,
    build_read_quality_comparison_rows,
    build_trimmed_quality_rows,
)
from cacao_inventory.schema import (  # noqa: E402
    CUTADAPT_SUMMARY_COLUMNS,
    RAW_READ_QUALITY_COLUMNS,
    READ_QUALITY_COMPARISON_COLUMNS,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--multiqc-fastqc", type=Path, required=True)
    parser.add_argument("--cutadapt-summary", type=Path, required=True)
    parser.add_argument("--raw-quality", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--comparison-output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    multiqc_rows = read_tsv(args.multiqc_fastqc, MULTIQC_FASTQC_COLUMNS)
    cutadapt_rows = read_tsv(args.cutadapt_summary, CUTADAPT_SUMMARY_COLUMNS)
    raw_rows = read_tsv(args.raw_quality, RAW_READ_QUALITY_COLUMNS)
    trimmed_rows = build_trimmed_quality_rows(multiqc_rows, cutadapt_rows)
    comparison_rows = build_read_quality_comparison_rows(raw_rows, trimmed_rows)
    write_tsv_atomic(args.output, trimmed_rows, RAW_READ_QUALITY_COLUMNS)
    write_tsv_atomic(
        args.comparison_output,
        comparison_rows,
        READ_QUALITY_COMPARISON_COLUMNS,
    )
    logging.info("Compared raw and trimmed FastQC for %d FASTQ files", len(trimmed_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
