#!/usr/bin/env python3
"""Separate organelles/known non-bacteria into auditable derived ASV tables.

Inputs: diversity config, independently validated DADA2 and primary taxonomy.
Outputs: data/processed candidate ASVs plus per-ASV/sample logs, hashes and SUCCESS.
Fails on altered upstream input, unsupported policy, mixed studies or empty samples.
Never modifies original DADA2/taxonomy tables or sequencing data.
"""

import argparse
import logging
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/dada2"))
sys.path.insert(0, str(ROOT / "scripts/taxonomy"))
from validate_outputs import validate_run
from validate_taxonomy import validate as validate_taxonomy
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, validate_policy, write_json
from cacao_inventory.bacterial_filter import DATA_PRODUCTS, REPORT_PRODUCTS, prepare_tables
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.provenance import file_records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/diversity.yaml")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    config = load_json_yaml(args.config)
    paths = validate_policy(config, ROOT)
    report = paths["filtering_dir"]
    state = begin_stage(report, ROOT)
    upstream = {"dada2": validate_run(paths["dada2_dir"], ROOT),
                "taxonomy": validate_taxonomy(paths["taxonomy_dir"], paths["dada2_dir"], ROOT)}
    tax_config = load_json_yaml(paths["taxonomy_dir"] / "config_snapshot.yaml")
    if config["filtering"]["minimum_bootstrap"] != tax_config["classification"]["primary_min_boot"] or config["filtering"]["tax_levels"] != tax_config["classification"]["tax_levels"]:
        raise ValueError("Filtering call differs from validated taxonomy policy")
    inputs = [args.config, *[paths["dada2_dir"] / name for name in ("asv_counts.tsv", "asv_sequences.tsv", "sample_metadata.tsv", "validation.json")],
              *[paths["taxonomy_dir"] / name for name in ("taxonomy.tsv", "provenance.json", "validation.json", "config_snapshot.yaml")],
              Path(__file__), ROOT / "scripts/taxonomy/validate_taxonomy.py", ROOT / "scripts/dada2/validate_outputs.py",
              *[ROOT / "python/cacao_inventory" / (name + ".py") for name in ("analysis_artifacts", "bacterial_filter", "config", "io", "provenance", "taxonomy_reference", "taxonomy_validation")],
              ROOT / "workflow/rules/pilot_diversity.smk", ROOT / "environment/conda-linux-64.lock"]
    before = file_records(inputs, ROOT)
    tables, summary = prepare_tables(paths["dada2_dir"], paths["taxonomy_dir"], config["filtering"])
    for name, rows in tables.items():
        target = paths["bacterial_dir"] if name in DATA_PRODUCTS else report
        write_tsv_atomic(target / name, rows, list(rows[0]))
    write_json(report / "summary.json", summary)
    products = [*[paths["bacterial_dir"] / name for name in DATA_PRODUCTS], *[report / name for name in REPORT_PRODUCTS]]
    finish_stage(report, ROOT, state, before, args.config, config, products, upstream,
                 {"Python": platform.python_version()}, sys.argv)
    logging.info("Retained %s ASVs and %s reads; %s ASVs excluded; all %s samples retained",
                 summary["retained_asvs"], summary["retained_reads"], summary["excluded_asvs"], summary["samples"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
