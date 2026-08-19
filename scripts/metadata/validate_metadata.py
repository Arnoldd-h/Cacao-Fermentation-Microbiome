#!/usr/bin/env python3
"""Validate versioned inventory tables without contacting external services."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_project_configuration  # noqa: E402
from cacao_inventory.io import read_inventory  # noqa: E402
from cacao_inventory.validation import validate_inventory  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "config.yaml")
    parser.add_argument("--datasets", type=Path, default=ROOT / "config" / "datasets.yaml")
    parser.add_argument("--input-dir", type=Path, default=ROOT / "metadata")
    parser.add_argument("--stamp", type=Path, help="Optional success stamp for Snakemake")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_config, _ = load_project_configuration(args.config, args.datasets)
    tables = read_inventory(args.input_dir)
    summary = validate_inventory(tables, project_config)
    if args.stamp:
        args.stamp.parent.mkdir(parents=True, exist_ok=True)
        args.stamp.write_text("metadata inventory valid\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
