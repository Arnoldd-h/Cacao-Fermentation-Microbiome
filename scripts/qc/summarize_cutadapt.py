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
from cacao_inventory.io import write_tsv_atomic  # noqa: E402
from cacao_inventory.pilot_workflow import load_pilot_manifest, processing_configuration
from cacao_inventory.provenance import qc_provenance
from cacao_inventory.qc import build_cutadapt_summary_rows  # noqa: E402
from cacao_inventory.schema import (  # noqa: E402
    CUTADAPT_SUMMARY_COLUMNS,
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
        help="Defaults to results/qc/pilot/cutadapt/<manifest study_id>",
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
    manifest = load_pilot_manifest(args.manifest)
    report_directory = args.report_directory or ROOT / "results" / "qc" / "pilot" / "cutadapt" / manifest[0]["study_id"]
    rows: list[dict[str, str]] = []
    inputs: list[Path] = [args.manifest]
    for sample in manifest:
        bioproject = sample["bioproject"]
        try:
            processing = processing_configuration(config, sample)
            primers = processing["primers"]
            parameters = processing["cutadapt"]
        except (KeyError, TypeError) as exc:
            raise ValueError(f"Missing Cutadapt configuration for {bioproject}") from exc
        run_accession = sample["run_accession"]
        report_path = report_directory / f"{run_accession}.cutadapt.json"
        inputs.append(report_path)
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
            inputs.append(ROOT / relative_path)
        inputs.extend(ROOT / report["input"][field] for field in ("path1", "path2"))
        rows.extend(
            build_cutadapt_summary_rows(
                report,
                sample,
                output_files,
                primers,
                parameters,
                root=ROOT,
                report_path=report_path,
            )
        )

    write_tsv_atomic(args.output, rows, CUTADAPT_SUMMARY_COLUMNS)
    qc_provenance(
        args.output, root=ROOT, inputs=inputs, tools=("cutadapt",),
        config=args.config, command=[sys.executable, *sys.argv],
    )
    logging.info("Recorded Cutadapt metrics for %d trimmed FASTQ files", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
