#!/usr/bin/env python3
"""Export the longitudinal hierarchy without modifying abundance data.

Inputs: validated studies.tsv/samples.tsv and config.analysis_design.
Outputs: results/tables/analysis_units.tsv, analysis_design.tsv and provenance.
Fails on missing study/batch/time, duplicate runs or conflicting inclusion.
No inferential model is fitted and no FASTQ are downloaded.
"""

from __future__ import annotations

import argparse
import logging
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_json_yaml  # noqa: E402
from cacao_inventory.design import DESIGN_COLUMNS, UNIT_COLUMNS, build_analysis_units  # noqa: E402
from cacao_inventory.io import read_inventory, write_tsv_atomic  # noqa: E402
from cacao_inventory.provenance import write_provenance  # noqa: E402
from cacao_inventory.validation import validate_inventory  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/config.yaml")
    parser.add_argument("--input-dir", type=Path, default=ROOT / "metadata")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/tables")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    config = load_json_yaml(args.config)
    tables = read_inventory(args.input_dir)
    validate_inventory(tables, config)
    units, summaries = build_analysis_units(tables["samples"], tables["studies"], config)
    unit_path = args.output_dir / "analysis_units.tsv"
    summary_path = args.output_dir / "analysis_design.tsv"
    write_tsv_atomic(unit_path, units, UNIT_COLUMNS)
    write_tsv_atomic(summary_path, summaries, DESIGN_COLUMNS)
    write_provenance(
        args.output_dir / "analysis_design.provenance.json", root=ROOT,
        inputs=[args.input_dir / "samples.tsv", args.input_dir / "studies.tsv",
                Path(__file__), ROOT / "python/cacao_inventory/design.py"],
        outputs=[unit_path, summary_path], config_paths=[args.config],
        software_versions={"Python": platform.python_version()}, command=sys.argv,
        parameters=config["analysis_design"],
    )
    logging.info("Mapped %d runs from %d studies into longitudinal units", len(units), len(summaries))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
