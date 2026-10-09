#!/usr/bin/env python3
"""Independently check bacterial filtering decisions and exact count conservation.

Inputs: configured source/derived ASVs, per-ASV/sample logs and provenance.
Output: JSON validation summary. Fails on unsupported exclusions, altered reads,
missing samples/ASVs, changed metadata, inconsistent summaries or checksum errors.
"""

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.analysis_artifacts import read_counts, read_table, unique_rows, validate_artifacts, validate_policy, write_json
from cacao_inventory.bacterial_filter import DATA_PRODUCTS, REPORT_PRODUCTS
from cacao_inventory.config import load_json_yaml


def validate_tables(config: dict, paths: dict) -> dict:
    source, features = read_counts(paths["dada2_dir"] / "asv_counts.tsv")
    derived, retained = read_counts(paths["bacterial_dir"] / "asv_counts.tsv")
    originals = unique_rows(source, ("study_id", "sample_id"))
    processed = unique_rows(derived, ("study_id", "sample_id"))
    if set(processed) != set(originals):
        raise ValueError("Sample exclusion or identity change")
    settings = config["filtering"]
    ranks = settings["tax_levels"]
    calls = unique_rows(read_table(paths["taxonomy_dir"] / "taxonomy.tsv", ["study_id", "asv_id", *ranks]), ("study_id", "asv_id"))
    log = unique_rows(read_table(paths["filtering_dir"] / "asv_filter_log.tsv", ["study_id", "asv_id", "retained", "reason", "total_reads", *ranks]), ("study_id", "asv_id"))
    study = source[0]["study_id"]
    if set(calls) != {(study, feature) for feature in features} or set(log) != set(calls):
        raise ValueError("ASV audit membership mismatch")
    expected_retained = set()
    for key, lineage in calls.items():
        organs = [label for label in settings["organelle_labels"] if any(lineage[rank].casefold() == label.casefold() for rank in ranks)]
        non_bacterial = bool(lineage["Kingdom"]) and lineage["Kingdom"] != settings["bacterial_kingdom"] and settings["exclude_known_non_bacterial"]
        reason = "organelle:" + ";".join(organs) if organs else "non_bacterial_kingdom:" + lineage["Kingdom"] if non_bacterial else "retained_bacterial" if lineage["Kingdom"] else "retained_unclassified_kingdom"
        keep = not organs and not non_bacterial
        if keep:
            expected_retained.add(key[1])
        if log[key]["retained"] != ("TRUE" if keep else "FALSE") or log[key]["reason"] != reason or any(log[key][rank] != lineage[rank] for rank in ranks):
            raise ValueError("Unsupported ASV exclusion or altered taxonomy")
        if log[key]["total_reads"] != str(sum(int(row[key[1]]) for row in source)):
            raise ValueError("Audit read total mismatch")
    if set(retained) != expected_retained:
        raise ValueError("Derived ASV membership mismatch")
    retention = unique_rows(read_table(paths["filtering_dir"] / "sample_retention.tsv", ["study_id", "sample_id", "input_reads", "retained_reads", "excluded_reads", "retained_fraction"]), ("study_id", "sample_id"))
    if set(retention) != set(originals):
        raise ValueError("Sample audit membership mismatch")
    for key, row in processed.items():
        if any(row[feature] != originals[key][feature] for feature in retained):
            raise ValueError("Retained count changed")
        total = sum(int(originals[key][feature]) for feature in features)
        kept = sum(int(row[feature]) for feature in retained)
        audit = retention[key]
        if not kept or any(audit[field] != str(value) for field, value in (("input_reads", total), ("retained_reads", kept), ("excluded_reads", total - kept))) or not math.isclose(float(audit["retained_fraction"]), kept / total, abs_tol=1e-10):
            raise ValueError("Sample retention balance mismatch")
    for name, directory, identifiers in (("asv_sequences.tsv", paths["dada2_dir"], ("study_id", "asv_id")),
                                        ("taxonomy.tsv", paths["taxonomy_dir"], ("study_id", "asv_id")),
                                        ("sample_metadata.tsv", paths["dada2_dir"], ("study_id", "sample_id"))):
        original = unique_rows(read_table(directory / name, list(identifiers)), identifiers)
        observed = unique_rows(read_table(paths["bacterial_dir"] / name, list(identifiers)), identifiers)
        expected = {key: value for key, value in original.items() if identifiers[1] == "sample_id" or key[1] in retained}
        if observed != expected:
            raise ValueError("Derived sequence/taxonomy/metadata changed")
    summary = {"study_id": study, "samples": len(source), "input_asvs": len(features), "retained_asvs": len(retained),
        "excluded_asvs": len(features) - len(retained), "input_reads": sum(sum(int(row[f]) for f in features) for row in source),
        "retained_reads": sum(sum(int(row[f]) for f in retained) for row in derived),
        "excluded_reads": sum(sum(int(row[f]) for f in features if f not in retained) for row in source),
        "unclassified_kingdom_retained_asvs": sum(not calls[study, f]["Kingdom"] for f in retained), "sample_exclusions": 0}
    if json.loads((paths["filtering_dir"] / "summary.json").read_text()) != summary:
        raise ValueError("Filtering summary mismatch")
    return summary


def validate(config_path: Path, root: Path) -> dict:
    config = load_json_yaml(config_path)
    paths = validate_policy(config, root)
    report = paths["filtering_dir"]
    expected = [*[paths["bacterial_dir"] / name for name in DATA_PRODUCTS], *[report / name for name in REPORT_PRODUCTS],
                report / "config_snapshot.yaml", report / "input_checksums.json"]
    summary = validate_artifacts(report, root, expected, config)
    summary.update(validate_tables(config, paths))
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
