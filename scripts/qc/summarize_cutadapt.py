#!/usr/bin/env python3
"""Build the validated per-read Cutadapt summary for the pilot.

Inputs: config/config.yaml, metadata/pilot_manifest.tsv, per-run Cutadapt JSON
reports, and paired trimmed FASTQ under data/interim. Output:
results/qc/pilot/cutadapt_summary.tsv by default. Fails when paired counts,
paths, parameters, or output files are inconsistent.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_json_yaml  # noqa: E402
from cacao_inventory.io import read_tsv, write_tsv_atomic  # noqa: E402
from cacao_inventory.qc import build_cutadapt_summary_rows  # noqa: E402
from cacao_inventory.schema import (  # noqa: E402
    CUTADAPT_SUMMARY_COLUMNS,
    PILOT_MANIFEST_COLUMNS,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "config.yaml")
    parser.add_argument(
        "--manifest", type=Path, default=ROOT / "metadata" / "pilot_manifest.tsv"
    )
    parser.add_argument(
        "--report-directory",
        type=Path,
        default=ROOT / "results" / "qc" / "pilot" / "cutadapt",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "qc" / "pilot" / "cutadapt_summary.tsv",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    config = load_json_yaml(args.config)
    manifest = read_tsv(args.manifest, PILOT_MANIFEST_COLUMNS)
    rows: list[dict[str, str]] = []
    for sample in manifest:
        bioproject = sample["bioproject"]
        try:
            processing = config["amplicon_processing"][bioproject]
            primers = processing["primers"]
            parameters = processing["cutadapt"]
        except (KeyError, TypeError) as exc:
            raise ValueError(f"Missing Cutadapt configuration for {bioproject}") from exc
        run_accession = sample["run_accession"]
        report_path = args.report_directory / f"{run_accession}.cutadapt.json"
        with report_path.open(encoding="utf-8") as handle:
            report = json.load(handle)
        relative_directory = (
            Path("data")
            / "interim"
            / sample["study_id"]
            / "pilot"
            / run_accession
        )
        output_files = {
            direction: (relative_directory / f"{run_accession}_{suffix}.fastq.gz").as_posix()
            for direction, suffix in (("R1", "1"), ("R2", "2"))
        }
        for relative_path in output_files.values():
            if not (ROOT / relative_path).is_file():
                raise FileNotFoundError(f"Missing trimmed FASTQ: {relative_path}")
        rows.extend(
            build_cutadapt_summary_rows(
                report,
                sample,
                output_files,
                primers,
                parameters,
            )
        )

    write_tsv_atomic(args.output, rows, CUTADAPT_SUMMARY_COLUMNS)
    logging.info("Recorded Cutadapt metrics for %d trimmed FASTQ files", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
