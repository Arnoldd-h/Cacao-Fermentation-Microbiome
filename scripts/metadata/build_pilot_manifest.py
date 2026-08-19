#!/usr/bin/env python3
"""Build a six-sample pilot manifest without downloading FASTQ.

Inputs: config/config.yaml, results/tables/pilot_dataset_selection.tsv and the
versioned runs/samples inventory. Output: metadata/pilot_manifest.tsv.
Fails if the selected study lacks enough independent batches or complete paired
FASTQ paths, byte counts and MD5 checksums in any required temporal stage.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_json_yaml  # noqa: E402
from cacao_inventory.io import read_inventory, read_tsv, write_tsv_atomic  # noqa: E402
from cacao_inventory.pilot import build_pilot_manifest  # noqa: E402
from cacao_inventory.schema import PILOT_MANIFEST_COLUMNS, PILOT_SELECTION_COLUMNS  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "config.yaml")
    parser.add_argument("--input-dir", type=Path, default=ROOT / "metadata")
    parser.add_argument(
        "--selection",
        type=Path,
        default=ROOT / "results" / "tables" / "pilot_dataset_selection.tsv",
    )
    parser.add_argument(
        "--output", type=Path, default=ROOT / "metadata" / "pilot_manifest.tsv"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_json_yaml(args.config)
    if "pilot_vertical_slice" not in config:
        raise ValueError("config.yaml must define pilot_vertical_slice")
    selection = read_tsv(args.selection, PILOT_SELECTION_COLUMNS)
    selected = [row for row in selection if row["selected_as_pilot"] == "true"]
    if len(selected) != 1:
        raise ValueError(f"Expected one selected pilot study, found {len(selected)}")
    inventory = read_inventory(args.input_dir)
    rows = build_pilot_manifest(
        selected[0]["study_id"],
        inventory["runs"],
        inventory["samples"],
        config["pilot_vertical_slice"],
    )
    write_tsv_atomic(args.output, rows, PILOT_MANIFEST_COLUMNS)
    print(
        json.dumps(
            {
                "estimated_fastq_bytes": sum(int(row["estimated_bytes_total"]) for row in rows),
                "pilot_bioproject": selected[0]["bioproject"],
                "samples": len(rows),
                "stages": sorted({row["fermentation_stage"] for row in rows}),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
