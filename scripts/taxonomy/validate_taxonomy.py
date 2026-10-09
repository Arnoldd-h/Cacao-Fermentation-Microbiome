#!/usr/bin/env python3
"""Independently validate completed taxonomy tables and original input hashes.

Inputs: taxonomy SUCCESS/provenance plus original DADA2 and reference inputs.
Output: JSON validation summary. Fails on altered inputs/artifacts, lost ASVs,
incorrect bootstrap masks, inconsistent flags or read/coverage totals.
"""

import argparse
import json
import logging
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.config import load_json_yaml
from cacao_inventory.provenance import file_records
from cacao_inventory.taxonomy_reference import repository_path, validate_taxonomy_config
from cacao_inventory.taxonomy_validation import validate_taxonomy_tables

REQUIRED = {"taxonomy.tsv", "taxonomy_unfiltered.tsv", "taxonomy_bootstraps.tsv", "taxonomy_sensitivity.tsv",
            "taxonomy_screening.tsv", "assignment_coverage.tsv", "sample_coverage.tsv", "config_snapshot.yaml",
            "software_versions.tsv", "session_info.txt", "warnings.txt", "input_checksums.json",
            "assignment_coverage.pdf", "assignment_coverage.svg", "assignment_coverage.png"}


def validate(directory: Path, input_directory: Path, root: Path) -> dict:
    directory, root = directory.resolve(), root.resolve()
    revision = (directory / "SUCCESS").read_text().strip()
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("Invalid taxonomy success marker")
    provenance = json.loads((directory / "provenance.json").read_text())
    if provenance["status"] != "success" or provenance["git_commit"] != revision:
        raise ValueError("Taxonomy success and provenance disagree")
    for section in ("inputs", "outputs"):
        records = provenance[section]
        paths = [repository_path(root, record["path"]) for record in records]
        if len(set(paths)) != len(paths):
            raise ValueError("Duplicate checksum path")
        if section == "outputs" and any(directory not in path.parents for path in paths):
            raise ValueError("Output path escapes taxonomy directory")
        if file_records(paths, root) != records:
            raise ValueError(f"Taxonomy {section} checksum mismatch")
    if not REQUIRED.issubset({Path(record["path"]).name for record in provenance["outputs"]}):
        raise ValueError("Provenance omits essential taxonomy artifacts")
    if json.loads((directory / "input_checksums.json").read_text()) != provenance["inputs"]:
        raise ValueError("Input checksum table disagrees with provenance")
    config = load_json_yaml(directory / "config_snapshot.yaml")
    validate_taxonomy_config(config, root)
    if provenance["classification"] != config["classification"] or provenance["reference"] != config["reference"]:
        raise ValueError("Executed taxonomy parameters disagree with configuration")
    summary = validate_taxonomy_tables(directory, input_directory, config)
    summary.update(status="valid", source_git_commit=revision,
                   artifacts_checked=len(provenance["outputs"]), inputs_checked=len(provenance["inputs"]))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=ROOT / "results/taxonomy/pilot")
    parser.add_argument("--input-dir", type=Path, default=ROOT / "results/dada2/pilot")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    result = validate(args.run_dir, args.input_dir, ROOT)
    document = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(".json.tmp")
        temporary.write_text(document, encoding="utf-8")
        temporary.replace(args.output)
    print(document, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
