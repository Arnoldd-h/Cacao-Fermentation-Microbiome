#!/usr/bin/env python3
"""Independently validate descriptive diversity formulas, PCA geometry and hashes.

Inputs: registered config, bacterial counts/metadata and R outputs/provenance.
Output: JSON validation summary. Fails on changed inputs, lost IDs, altered alpha,
incorrect CLR/distances, nonfinite values or PCA centering/variance/geometry errors.
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from run_diversity import PRODUCTS
from cacao_inventory.analysis_artifacts import validate_artifacts, validate_policy, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.diversity_validation import validate_diversity_tables


def validate(config_path: Path, root: Path) -> dict:
    config = load_json_yaml(config_path)
    paths = validate_policy(config, root)
    report = paths["diversity_dir"]
    expected = [*[report / name for name in PRODUCTS], report / "config_snapshot.yaml", report / "input_checksums.json"]
    summary = validate_artifacts(report, root, expected, config)
    summary.update(validate_diversity_tables(report, paths["bacterial_dir"], config))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/diversity.yaml")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    summary = validate(args.config, ROOT)
    if args.output:
        write_json(args.output, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
