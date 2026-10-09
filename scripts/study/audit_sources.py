#!/usr/bin/env python3
"""Verify full-study NCBI/ENA metadata before auditing FASTQ packaging.

Inputs: registered study scope, candidate rules and versioned runs.
Output: small source-audit JSON; no FASTQ downloads or inventory edits.
Fails on changed accession membership or download/design metadata.
"""
import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.analysis_artifacts import write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import read_inventory
from cacao_inventory.provenance import write_provenance
from cacao_inventory.sources import MetadataClient
from cacao_inventory.study_scope import compare_live_runs, validate_scope

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/full_study.yaml")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    cfg = load_json_yaml(args.config)
    validate_scope(cfg, ROOT)
    project = load_json_yaml(ROOT / cfg["project_config"])
    candidate = next(row for row in load_json_yaml(ROOT / cfg["datasets_config"])["candidates"] if row["bioproject"] == cfg["bioproject"])
    client = MetadataClient(project["inventory"])
    summary = client.fetch_bioproject_summaries([cfg["bioproject"]])[cfg["bioproject"]]
    audit = compare_live_runs(cfg, candidate, read_inventory(ROOT / cfg["inventory_dir"])["runs"], client.fetch_ena_runs(cfg["bioproject"]))
    audit["ncbi_bioproject"] = summary["project_acc"]
    audit["ncbi_title"] = summary.get("project_title", "")
    output = ROOT / cfg["plan_dir"] / "packaging_source_audit.json"
    write_json(output, audit)
    write_provenance(output.with_suffix(".provenance.json"), root=ROOT,
                     inputs=[ROOT / cfg["inventory_dir"] / "runs.tsv", ROOT / cfg["datasets_config"], Path(__file__), ROOT / "python/cacao_inventory/study_scope.py", ROOT / "python/cacao_inventory/sources.py"],
                     outputs=[output], config_paths=[args.config, ROOT / cfg["project_config"]],
                     software_versions={"Python": sys.version.split()[0]}, command=sys.argv)
    logging.info("NCBI/ENA audit passed: %d eligible runs", audit["eligible_runs"])

if __name__ == "__main__":
    main()
