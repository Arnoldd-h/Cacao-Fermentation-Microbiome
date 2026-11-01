#!/usr/bin/env python3
"""Fit registered within-study CLR trajectories with independent-batch inference.

Inputs: longitudinal config and independently validated bacterial counts/metadata.
Outputs: membership, OLS coefficients, global exact signs/BY tests, heterogeneity,
leave-one-batch-out, figures, independent R checks, hashes and SUCCESS/validation.
Fails on invalid/stale upstream data, partial designs, numerical disagreement,
changed inputs or R errors. No raw/previous-stage tables are modified.
--validate-only checks existing results without rewriting tables or rerunning R.
"""

import argparse
import json
import logging
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.provenance import file_records
from longitudinal_model import build_tables, table_columns
from validate_longitudinal import PRODUCTS, load_sources, model_summary, validate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/longitudinal.yaml")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    config = load_json_yaml(args.config)
    source, report, paths, upstream, counts, features, metadata = load_sources(config, ROOT)
    project = load_json_yaml(ROOT / "config/config.yaml")
    if config["method"]["random_seed"] != project["project"]["default_random_seed"]:
        raise ValueError("Longitudinal and project seeds disagree")
    if not args.validate_only:
        state = begin_stage(report, ROOT)
        inputs = [args.config, source, ROOT / "config/config.yaml", ROOT / "environment/conda-linux-64.lock",
                  ROOT / "workflow/longitudinal.Snakefile", *[path for path in Path(__file__).parent.glob("*") if path.suffix in (".py", ".R")],
                  *[paths["bacterial_dir"] / name for name in ("asv_counts.tsv", "sample_metadata.tsv", "taxonomy.tsv")],
                  *[paths["filtering_dir"] / name for name in ("SUCCESS", "validation.json", "provenance.json")],
                  *[ROOT / "scripts" / name for name in ("filtering/validate_bacterial_table.py", "taxonomy/validate_taxonomy.py", "taxonomy/seeded_policy.py")],
                  *[ROOT / "python/cacao_inventory" / (name + ".py") for name in ("analysis_artifacts", "bacterial_filter", "config", "io", "provenance", "taxonomy_reference", "taxonomy_validation")]]
        before = file_records(inputs, ROOT)
        logging.info("Fitting %s registered scenarios with %s independent batches", len(config["scenarios"]), config["cohort"]["expected_batches"])
        tables = build_tables(counts, features, metadata, config)
        for name, columns in table_columns(features).items():
            write_tsv_atomic(report / name, tables[name], columns)
        write_json(report / "summary.json", model_summary(config, features, tables))
        commands = [
            ["Rscript", "scripts/longitudinal/validate_model.R", "--config", str(args.config.resolve()), "--counts", str(paths["bacterial_dir"] / "asv_counts.tsv"),
             "--metadata", str(paths["bacterial_dir"] / "sample_metadata.tsv"), "--output-dir", str(report)],
            ["Rscript", "scripts/longitudinal/plot_trajectories.R", "--config", str(args.config.resolve()), "--input", str(report / "trajectory_plot_data.tsv"), "--output-dir", str(report)]]
        for command in commands:
            subprocess.run(command, cwd=ROOT, check=True)
        finish_stage(report, ROOT, state, before, args.config, config, [report / name for name in PRODUCTS], upstream,
                     {"Python": platform.python_version()}, commands)
    result = validate(args.config)
    write_json(report / "validation.json", result)
    logging.info("Longitudinal analysis validated: %s samples, %s observations, %s independent batches",
                 result["primary_samples"], result["primary_observations"], result["independent_batches"])
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
