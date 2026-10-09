#!/usr/bin/env python3
"""Describe validated candidate ASVs in R and record reproducible provenance.

Inputs: registered config and independently validated bacterial filtering stage.
Outputs: diversity tables/figures, snapshots/hashes, provenance and SUCCESS.
Fails on altered upstream data, unsupported methods, seed mismatch or R failure.
No inferential tests, sample aggregation, rarefaction or cross-study pooling.
"""

import argparse
import logging
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/filtering"))
from validate_bacterial_table import validate as validate_filtering
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, validate_policy
from cacao_inventory.bacterial_filter import DATA_PRODUCTS
from cacao_inventory.config import load_json_yaml
from cacao_inventory.provenance import file_records

PRODUCTS = ["alpha_diversity.tsv", "clr_coordinates.tsv", "aitchison_distances.tsv", "bray_curtis_distances.tsv",
            "pca_scores.tsv", "pca_variance.tsv", "software_versions.tsv", "session_info.txt", "warnings.txt",
            *[name + "." + extension for name in ("alpha_diversity", "aitchison_pca") for extension in ("pdf", "svg", "png")]]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/diversity.yaml")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    config = load_json_yaml(args.config)
    paths = validate_policy(config, ROOT)
    report = paths["diversity_dir"]
    state = begin_stage(report, ROOT)
    project = load_json_yaml(ROOT / "config/config.yaml")
    if config["analysis"]["random_seed"] != project["project"]["default_random_seed"] or project["analysis_design"]["pilot_inference"]:
        raise ValueError("Seed or pilot inference policy mismatch")
    upstream = validate_filtering(args.config, ROOT)
    inputs = [args.config, ROOT / "config/config.yaml", *[paths["bacterial_dir"] / name for name in DATA_PRODUCTS],
              *[paths["filtering_dir"] / name for name in ("validation.json", "provenance.json")], Path(__file__),
              ROOT / "scripts/diversity/describe_pilot.R", ROOT / "scripts/dada2/helpers.R", ROOT / "scripts/filtering/validate_bacterial_table.py",
              *[ROOT / "python/cacao_inventory" / (name + ".py") for name in ("analysis_artifacts", "bacterial_filter", "config", "io", "provenance", "taxonomy_reference")],
              ROOT / "workflow/rules/pilot_diversity.smk", ROOT / "environment/conda-linux-64.lock"]
    before = file_records(inputs, ROOT)
    command = ["Rscript", "scripts/diversity/describe_pilot.R", "--config", str(args.config.resolve())]
    subprocess.run(command, cwd=ROOT, check=True)
    finish_stage(report, ROOT, state, before, args.config, config, [report / name for name in PRODUCTS],
                 upstream, {"Python": platform.python_version()}, command)
    logging.info("Descriptive diversity completed for %s samples", upstream["samples"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
