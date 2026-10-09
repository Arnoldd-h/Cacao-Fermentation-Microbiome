#!/usr/bin/env python3
"""Audit exceptional ENA FASTQ packaging after the separate public-source audit.

Inputs: registered study, valid source audit and original inventory FASTQ records.
Outputs: small byte/MD5/gzip/read/header summaries of nonstandard run streams.
Downloads only those exceptional streams into immutable data/raw; never invents mates.
Fails on unaudited metadata, altered files, malformed FASTQ or network errors.
"""
import argparse
import gzip
import json
import logging
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.analysis_artifacts import write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.download import download_file, expand_pilot_manifest
from cacao_inventory.fastq import validate_fastq_file
from cacao_inventory.io import read_inventory
from cacao_inventory.provenance import file_records, write_provenance
from cacao_inventory.study_scope import validate_scope

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/full_study.yaml")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    cfg = load_json_yaml(args.config)
    validate_scope(cfg, ROOT)
    audit_path = ROOT / cfg["plan_dir"] / "packaging_source_audit.json"
    audit = json.loads(audit_path.read_text())
    provenance = json.loads(audit_path.with_suffix(".provenance.json").read_text())
    for section in ("config", "inputs", "outputs"):
        recorded = provenance[section]
        if file_records([ROOT / row["path"] for row in recorded], ROOT) != recorded:
            raise ValueError("Source audit checksum changed: " + section)
    runs = [row for row in read_inventory(ROOT / cfg["inventory_dir"])["runs"] if row["study_id"] == cfg["study_id"] and row["analysis_include"] == "true"]
    if audit["status"] != "valid" or set(audit["verified_run_accessions"]) != {row["run_accession"] for row in runs}:
        raise ValueError("Valid source audit is required before FASTQ inspection")
    rows = []
    inputs = [args.config, audit_path, audit_path.with_suffix(".provenance.json"), ROOT / cfg["inventory_dir"] / "runs.tsv", Path(__file__)]
    for run in runs:
        if len(run["fastq_ftp"].split(";")) == 2:
            continue
        row = {**run, "sample_id": run["sample_alias"], "fastq_ftp": ";".join("https://" + url if not url.startswith(("http://", "https://", "ftp://")) else url for url in run["fastq_ftp"].split(";"))}
        for record in expand_pilot_manifest([row]):
            path = ROOT / "data/raw" / record["target_relative"]
            download_file(record["source_url"], path, int(record["expected_bytes"]), record["expected_md5"], 1048576)
            metrics = validate_fastq_file(path)
            identifiers = Counter()
            mate_tokens = Counter()
            examples = []
            with gzip.open(path, "rt") as handle:
                while header := handle.readline():
                    if len(examples) < 4:
                        examples.append(header.strip())
                    fields = header.strip().split()
                    identifiers[fields[0]] += 1
                    suffix = fields[1].rsplit("/", 1)[-1] if len(fields) > 1 and "/" in fields[1] else "unknown"
                    mate_tokens[suffix if suffix in ("1", "2") else "unknown"] += 1
                    for _ in range(3):
                        handle.readline()
            rows.append({**record, **metrics, "distinct_header_ids": len(identifiers), "repeated_header_ids": sum(value > 1 for value in identifiers.values()),
                         "header_mate_suffixes": dict(mate_tokens), "header_examples": examples})
            inputs.append(path)
            logging.info("%s: %d reads, %d distinct header IDs", path.name, metrics["read_count"], len(identifiers))
    output = ROOT / cfg["plan_dir"] / "packaging_audit.json"
    write_json(output, {"status": "valid", "streams": rows})
    write_provenance(output.with_suffix(".provenance.json"), root=ROOT, inputs=inputs, outputs=[output],
                     config_paths=[args.config], software_versions={"Python": sys.version.split()[0]}, command=sys.argv)

if __name__ == "__main__":
    main()
