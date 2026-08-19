#!/usr/bin/env python3
"""Select the pilot dataset reproducibly from versioned inventory metadata.

Inputs: config/config.yaml and metadata/studies.tsv, runs.tsv and samples.tsv.
Output: results/tables/pilot_dataset_selection.tsv.
Fails on missing schemas, invalid selection settings or no eligible primary study.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_json_yaml  # noqa: E402
from cacao_inventory.io import read_inventory, write_tsv_atomic  # noqa: E402
from cacao_inventory.schema import PILOT_SELECTION_COLUMNS  # noqa: E402
from cacao_inventory.selection import build_pilot_selection  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "config.yaml")
    parser.add_argument("--input-dir", type=Path, default=ROOT / "metadata")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "tables" / "pilot_dataset_selection.tsv",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_json_yaml(args.config)
    if "pilot_selection" not in config:
        raise ValueError("config.yaml must define pilot_selection")
    inventory = read_inventory(args.input_dir)
    rows = build_pilot_selection(
        inventory["studies"], inventory["runs"], inventory["samples"], config["pilot_selection"]
    )
    selected = [row for row in rows if row["selected_as_pilot"] == "true"]
    if len(selected) != 1:
        raise ValueError(f"Expected exactly one selected pilot, found {len(selected)}")
    write_tsv_atomic(args.output, rows, PILOT_SELECTION_COLUMNS)
    print(
        json.dumps(
            {
                "eligible_studies": sum(row["eligible_primary"] == "true" for row in rows),
                "pilot_bioproject": selected[0]["bioproject"],
                "pilot_study_id": selected[0]["study_id"],
                "studies_compared": len(rows),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
