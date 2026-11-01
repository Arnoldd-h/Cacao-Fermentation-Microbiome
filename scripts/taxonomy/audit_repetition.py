#!/usr/bin/env python3
"""Compare complete taxonomy tables from two chronological classifications.

Inputs: two completed taxonomy directories, including SUCCESS and config snapshots.
Output: JSON SHA-256 and cell-level differences for seven scientific tables.
Fails on missing files, duplicate/changed row keys or (with --require-identical)
any byte difference. Tables and classifications are never edited or selected.
"""

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

TABLES = ("taxonomy.tsv", "taxonomy_unfiltered.tsv", "taxonomy_bootstraps.tsv",
          "taxonomy_sensitivity.tsv", "taxonomy_screening.tsv", "assignment_coverage.tsv", "sample_coverage.tsv")


def read_table(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames
        keys = [column for column in ("study_id", "asv_id", "sample_id", "min_boot", "rank") if column in columns]
        rows = {}
        for row in reader:
            key = tuple(row[column] for column in keys)
            if key in rows or None in row or any(value is None for value in row.values()):
                raise ValueError(f"Duplicate row or invalid schema: {path}")
            rows[key] = row
    return columns, keys, rows


def compare_tables(baseline: Path, repeat: Path) -> dict:
    tables = {}
    for name in TABLES:
        left, right = baseline / name, repeat / name
        columns, keys, before = read_table(left)
        other_columns, other_keys, after = read_table(right)
        if columns != other_columns or keys != other_keys or before.keys() != after.keys():
            raise ValueError(f"Classification row keys/schema changed: {name}")
        differences = []
        changed_rows = set()
        for key, row in before.items():
            for column in columns:
                if row[column] != after[key][column]:
                    changed_rows.add(key)
                    differences.append({"key": dict(zip(keys, key)), "column": column,
                                        "before": row[column], "after": after[key][column]})
        before_hash = hashlib.sha256(left.read_bytes()).hexdigest()
        after_hash = hashlib.sha256(right.read_bytes()).hexdigest()
        tables[name] = {"baseline_sha256": before_hash, "repeat_sha256": after_hash,
                        "identical": before_hash == after_hash, "rows": len(before),
                        "changed_rows": len(changed_rows), "changed_cells": len(differences),
                        "differences": differences}
    return tables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--repeat", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-identical", action="store_true")
    args = parser.parse_args()
    if args.baseline.resolve() == args.repeat.resolve():
        raise ValueError("Independent execution directories are required")
    tables = compare_tables(args.baseline, args.repeat)
    record = {"compared_at_utc": datetime.now(timezone.utc).isoformat(), "tables": tables,
              "all_identical": all(row["identical"] for row in tables.values()),
              "executions": {}}
    for label, directory in (("baseline", args.baseline), ("repeat", args.repeat)):
        provenance_path = directory / "provenance.json"
        provenance = json.loads(provenance_path.read_text())
        record["executions"][label] = {"directory": str(directory),
            "git_commit": (directory / "SUCCESS").read_text().strip(),
            "effective_threads": provenance.get("effective_threads"),
            "started_at_utc": provenance.get("started_at_utc"),
            "finished_at_utc": provenance.get("finished_at_utc"),
            "provenance_sha256": hashlib.sha256(provenance_path.read_bytes()).hexdigest(),
            "config": json.loads((directory / "config_snapshot.yaml").read_text())}
    if args.require_identical:
        inputs = [sorted(json.loads((directory / "input_checksums.json").read_text()), key=lambda row: row["path"])
                  for directory in (args.baseline, args.repeat)]
        if inputs[0] != inputs[1] or record["executions"]["baseline"]["config"] != record["executions"]["repeat"]["config"]:
            raise ValueError("Repetition inputs/configuration differ")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"all_identical": record["all_identical"], "tables": {
        name: {key: row[key] for key in ("identical", "changed_rows", "changed_cells")} for name, row in tables.items()}}, indent=2))
    if args.require_identical and not record["all_identical"]:
        raise SystemExit("Taxonomy repetition differs; see complete audit")


if __name__ == "__main__":
    main()
