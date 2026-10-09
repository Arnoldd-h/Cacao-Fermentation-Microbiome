#!/usr/bin/env python3
"""Trim a manifest run with its own BioProject primer configuration.

Inputs: paired data/raw FASTQ, pilot manifest and processing configuration.
Outputs: paired data/interim FASTQ and Cutadapt JSON. Fails on unsupported
methodological flags, inconsistent identifiers, Cutadapt failure or report mismatch.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_json_yaml
from cacao_inventory.pilot_workflow import fastq_path, load_pilot_manifest, manifest_row, processing_configuration
from cacao_inventory.qc import build_cutadapt_summary_rows, cutadapt_command


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--run-accession", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    row = manifest_row(load_pilot_manifest(args.manifest), args.run_accession, args.study_id)
    processing = processing_configuration(load_json_yaml(args.config), row)
    outputs = {direction: fastq_path(row, suffix, "trimmed") for direction, suffix in (("R1", "1"), ("R2", "2"))}
    command = cutadapt_command(row, processing["primers"], processing["cutadapt"], outputs, str(args.report), args.threads)
    for path in [*(ROOT / value for value in outputs.values()), args.report]:
        path.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Trimming %s using %s configuration", args.run_accession, row["bioproject"])
    subprocess.run(command, cwd=ROOT, check=True)
    report = json.loads(args.report.read_text(encoding="utf-8"))
    build_cutadapt_summary_rows(report, row, outputs, processing["primers"], processing["cutadapt"], root=ROOT, report_path=args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
