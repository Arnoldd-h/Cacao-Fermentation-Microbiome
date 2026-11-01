#!/usr/bin/env python3
"""Validate complete within-study longitudinal outputs without rewriting tables.

Inputs: registered config, validated bacterial sources, result tables/provenance.
Output: validation.json (CLI) and a semantic/hash summary. Fails on stale inputs,
altered membership/coefficients/probabilities, R-check failures or invalid figures.
"""

import argparse
import json
import logging
import math
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/filtering"))
from validate_bacterial_table import validate as validate_filtering
from cacao_inventory.analysis_artifacts import read_counts, read_table, validate_artifacts, validate_policy, write_json
from cacao_inventory.config import load_json_yaml
from longitudinal_model import IDENTITY, build_tables, table_columns
from longitudinal_policy import validate_longitudinal_policy

TABLE_NAMES = ["observations.tsv", "sample_membership.tsv", "batch_models.tsv", "slopes.tsv", "intercepts.tsv", "global_tests.tsv",
               "direction_cosines.tsv", "leave_one_batch_out.tsv", "sign_statistics.tsv", "trajectory_plot_data.tsv"]
PRODUCTS = [*TABLE_NAMES, "summary.json", "independent_validation.tsv", "software_versions.tsv", "session_info.txt",
            *["community_trajectories." + extension for extension in ("pdf", "svg", "png")]]


def load_sources(config, root):
    source, report = validate_longitudinal_policy(config, root)
    paths = validate_policy(load_json_yaml(source), root)
    if any(report == path or report in path.parents or path in report.parents for path in paths.values()):
        raise ValueError("Longitudinal outputs overlap upstream stage directories")
    upstream = validate_filtering(source, root)
    counts, features = read_counts(paths["bacterial_dir"] / "asv_counts.tsv")
    metadata = read_table(paths["bacterial_dir"] / "sample_metadata.tsv", IDENTITY)
    return source, report, paths, upstream, counts, features, metadata


def model_summary(config, features, tables):
    primary = next(row for row in tables["global_tests.tsv"] if row["scenario_id"] == config["primary_scenario"])
    return {"study_id": config["cohort"]["study_id"], "scope": config["method"]["scope"],
            "source_samples": config["cohort"]["expected_samples"], "independent_batches": config["cohort"]["expected_batches"],
            "asvs": len(features), "primary_scenario": config["primary_scenario"], "primary_observations": primary["n_observations"],
            "primary_samples": primary["n_samples"], "primary_samples_not_used": config["cohort"]["expected_samples"] - primary["n_samples"],
            "scenarios": len(config["scenarios"]), "sign_patterns_per_scenario": primary["sign_patterns"],
            "minimum_attainable_p": primary["minimum_attainable_p"], "fdr_method": config["method"]["fdr_method"],
            "null_assumptions": config["method"]["null_assumptions"], "cross_study_inference": False, "pilot_inference": False,
            "random_sampling": False, "random_seed": config["method"]["random_seed"]}


def compare_table(path, expected, columns, tolerance):
    observed = read_table(path, columns)
    if len(observed) != len(expected):
        raise ValueError(f"Longitudinal row count mismatch: {path.name}")
    for row, target in zip(observed, expected):
        if set(row) != set(target):
            raise ValueError(f"Longitudinal columns mismatch: {path.name}")
        for field, value in target.items():
            if isinstance(value, (int, float)):
                try:
                    number = float(row[field])
                except ValueError as error:
                    raise ValueError(f"Invalid numeric output {path.name}/{field}") from error
                if not math.isfinite(number) or not math.isclose(number, value, rel_tol=tolerance, abs_tol=tolerance):
                    raise ValueError(f"Longitudinal formula mismatch: {path.name}/{field}")
            elif row[field] != value:
                raise ValueError(f"Longitudinal identity/policy mismatch: {path.name}/{field}")


def validate(config_path, root=ROOT):
    config = load_json_yaml(config_path)
    _, report, _, _, counts, features, metadata = load_sources(config, root)
    record = validate_artifacts(report, root, [*[report / name for name in PRODUCTS], report / "config_snapshot.yaml", report / "input_checksums.json"], config)
    tables = build_tables(counts, features, metadata, config)
    for name, columns in table_columns(features).items():
        compare_table(report / name, tables[name], columns, config["method"]["validation_tolerance"])
    summary = model_summary(config, features, tables)
    if json.loads((report / "summary.json").read_text()) != summary:
        raise ValueError("Longitudinal summary disagrees with validated data")
    checks = read_table(report / "independent_validation.tsv", ["scenario_id", "status", "asv_coefficients_checked", "observations_checked", "batches_checked"])
    if len(checks) != len(config["scenarios"]) or [row["scenario_id"] for row in checks] != [row["id"] for row in config["scenarios"]]:
        raise ValueError("Independent R validation scenario membership mismatch")
    for check, global_row in zip(checks, tables["global_tests.tsv"]):
        if check["status"] != "passed" or int(check["asv_coefficients_checked"]) != 2 * len(features) * global_row["n_batches"] or int(check["observations_checked"]) != global_row["n_observations"] or int(check["batches_checked"]) != global_row["n_batches"]:
            raise ValueError("Independent R validation incomplete or failed")
    with (report / "community_trajectories.png").open("rb") as handle:
        header = handle.read(24)
    figures = config["figures"]
    dimensions = tuple(round(figures[key] * figures["png_dpi"]) for key in ("width_inches", "height_inches"))
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or struct.unpack(">II", header[16:24]) != dimensions:
        raise ValueError("Longitudinal PNG dimensions invalid")
    if not (report / "community_trajectories.pdf").read_bytes().startswith(b"%PDF-") or "<svg" not in (report / "community_trajectories.svg").read_text():
        raise ValueError("Invalid longitudinal vector figures")
    return {**record, **summary, "independent_r_scenarios_checked": len(checks)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/longitudinal.yaml")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    result = validate(args.config)
    _, report = validate_longitudinal_policy(load_json_yaml(args.config), ROOT)
    write_json(report / "validation.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
