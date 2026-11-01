#!/usr/bin/env python3
"""Plot validated alpha diversity against observed hours within batch/stratum.

Inputs: temporal reporting profile and validated diversity alpha/metadata tables.
Outputs: temporal_diversity.pdf/.svg/.png (300 dpi), exact plot_data.tsv,
session info, summary, configuration, SHA-256 provenance and validation.json.
Fails on stale upstream inputs, mixed studies, missing/duplicate identities or
batch/hour/stratum observations, invalid numbers, unsupported policy or R failure.
--validate-only checks existing outputs without rewriting them.
No sample aggregation, interpolation, model fitting or inferential tests.
"""

import argparse
import json
import logging
import math
import platform
import re
import struct
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/diversity"))
from validate_diversity import validate as validate_source
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, read_table, unique_rows, validate_artifacts, validate_policy, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.provenance import file_records
from cacao_inventory.taxonomy_reference import repository_path

IDENTITY = ["study_id", "sample_id", "run_accession", "fermentation_batch", "fermentation_hours", "relative_time", "fermentation_stage", "sampling_stratum"]
METRICS = ["observed_asvs", "shannon_effective", "library_reads"]
COLUMNS = [*IDENTITY, *METRICS]
PRODUCTS = ["plot_data.tsv", "summary.json", "session_info.txt", *["temporal_diversity." + extension for extension in ("pdf", "svg", "png")]]


def temporal_profile(config, analysis, root):
    expected = {"x_field": "fermentation_hours", "group_fields": ["fermentation_batch", "sampling_stratum"],
                "sample_aggregation": "none", "inference": False, "break_at_missing_batch_times": True}
    if any(config.get(key) != value or (isinstance(value, bool) and type(config.get(key)) is not bool) for key, value in expected.items()):
        raise ValueError("Unsupported temporal plotting policy")
    if [row["field"] for row in config["metrics"]] != METRICS or any(not row["label"] for row in config["metrics"]):
        raise ValueError("Unexpected temporal metrics or missing labels")
    profile = config["analyses"].get(analysis)
    if not profile or any(type(profile[key]) is not int or profile[key] < 1 for key in ("expected_samples", "expected_batches")) or type(profile["connect_points"]) is not bool or not profile["scope_label"]:
        raise ValueError("Invalid temporal profile or scope")
    figures = config["figures"]
    if type(figures["png_dpi"]) is not int or figures["png_dpi"] < 300:
        raise ValueError("Temporal PNG requires at least 300 dpi")
    for key in ("width_inches", "height_inches"):
        value = figures[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError("Invalid temporal figure dimensions")
    if not figures["stratum_colors"] or any(not re.fullmatch(r"#[0-9a-fA-F]{6}", color) for color in figures["stratum_colors"]):
        raise ValueError("Invalid stratum colors")
    output = repository_path(root, profile["output_dir"])
    if not output.is_relative_to(root.resolve() / "results/temporal"):
        raise ValueError("Temporal outputs must remain under results/temporal")
    return profile, repository_path(root, profile["diversity_config"]), output


def prepare_plot_rows(rows, expected_samples, expected_batches):
    samples = unique_rows(rows, ("study_id", "sample_id"))
    if len(samples) != expected_samples or len({row["study_id"] for row in rows}) != 1 or len({row["fermentation_batch"] for row in rows}) != expected_batches:
        raise ValueError("Temporal sample/study/batch scope mismatch")
    if len({row["run_accession"] for row in rows}) != len(rows):
        raise ValueError("Run accession repeated across temporal samples")
    observations = set()
    for row in rows:
        if any(not row[field] for field in IDENTITY):
            raise ValueError("Missing temporal identity or metadata; no imputation")
        hours = float(row["fermentation_hours"])
        if not math.isfinite(hours) or hours < 0:
            raise ValueError("Invalid fermentation hours")
        key = (row["fermentation_batch"], hours, row["sampling_stratum"])
        if key in observations:
            raise ValueError("Duplicate batch/hour/stratum; aggregation is not registered")
        observations.add(key)
        for field in METRICS:
            value = float(row[field])
            if not math.isfinite(value) or value <= 0:
                raise ValueError("Invalid temporal metric")
        if any(not row[field].isdigit() for field in ("library_reads", "observed_asvs")):
            raise ValueError("Depth and ASV richness must be integer counts")
    ordered = sorted(rows, key=lambda row: (row["fermentation_batch"], float(row["fermentation_hours"]), row["sampling_stratum"], row["sample_id"]))
    return [{field: row[field] for field in COLUMNS} for row in ordered]


def table_summary(rows, connect_points):
    return {"study_id": rows[0]["study_id"], "samples": len(rows),
            "batches": len({row["fermentation_batch"] for row in rows}),
            "series": len({(row["fermentation_batch"], row["sampling_stratum"]) for row in rows}),
            "metrics": METRICS, "connect_points": connect_points,
            "sample_aggregation": "none", "inference": False}


def validate_figures(report, figures):
    if not (report / "temporal_diversity.pdf").read_bytes().startswith(b"%PDF-"):
        raise ValueError("Invalid temporal PDF")
    if "<svg" not in (report / "temporal_diversity.svg").read_text():
        raise ValueError("Invalid temporal SVG")
    with (report / "temporal_diversity.png").open("rb") as handle:
        header = handle.read(24)
    expected = tuple(round(figures[key] * figures["png_dpi"]) for key in ("width_inches", "height_inches"))
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or struct.unpack(">II", header[16:24]) != expected:
        raise ValueError("Temporal PNG dimensions/DPI contract mismatch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/temporal_reporting.yaml")
    parser.add_argument("--analysis", required=True, help="Analysis key in the temporal reporting configuration")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    config = load_json_yaml(args.config)
    profile, source_config, report = temporal_profile(config, args.analysis, ROOT)
    paths = validate_policy(load_json_yaml(source_config), ROOT)
    if any(report == path or report in path.parents or path in report.parents for path in paths.values()):
        raise ValueError("Temporal output must not contain upstream stage directories")
    upstream = validate_source(source_config, ROOT)
    source = paths["diversity_dir"]
    rows = prepare_plot_rows(read_table(source / "alpha_diversity.tsv", COLUMNS), profile["expected_samples"], profile["expected_batches"])
    if len({row["sampling_stratum"] for row in rows}) > len(config["figures"]["stratum_colors"]):
        raise ValueError("Insufficient colors for observed strata")
    summary = table_summary(rows, profile["connect_points"])
    if not args.validate_only:
        state = begin_stage(report, ROOT)
        inputs = [args.config, source_config, source / "alpha_diversity.tsv", source / "validation.json", source / "provenance.json",
                  Path(__file__), ROOT / "scripts/reporting/plot_temporal_diversity.R", ROOT / "scripts/diversity/validate_diversity.py",
                  ROOT / "scripts/diversity/run_diversity.py", ROOT / "python/cacao_inventory/diversity_validation.py",
                  *[ROOT / "python/cacao_inventory" / (name + ".py") for name in ("analysis_artifacts", "config", "io", "provenance", "taxonomy_reference")],
                  ROOT / "environment/conda-linux-64.lock"]
        before = file_records(inputs, ROOT)
        write_tsv_atomic(report / "plot_data.tsv", rows, COLUMNS)
        write_json(report / "summary.json", summary)
        command = ["Rscript", "scripts/reporting/plot_temporal_diversity.R", "--config", str(args.config.resolve()),
                   "--analysis", args.analysis, "--input", str(report / "plot_data.tsv"), "--output-dir", str(report)]
        subprocess.run(command, cwd=ROOT, check=True)
        finish_stage(report, ROOT, state, before, args.config, config, [report / name for name in PRODUCTS], upstream,
                     {"Python": platform.python_version()}, command)
    result = validate_artifacts(report, ROOT, [*[report / name for name in PRODUCTS], report / "config_snapshot.yaml", report / "input_checksums.json"], config)
    if read_table(report / "plot_data.tsv", COLUMNS) != rows or json.loads((report / "summary.json").read_text()) != summary:
        raise ValueError("Temporal plot data or summary changed")
    validate_figures(report, config["figures"])
    result.update(summary)
    if not args.validate_only:
        write_json(report / "validation.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
