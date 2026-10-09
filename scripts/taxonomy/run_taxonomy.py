#!/usr/bin/env python3
"""Validate upstream DADA2 and pinned reference, then classify one study in R.

Inputs: taxonomy config, validated DADA2 directory and reference with SHA-256 sidecar.
Outputs: taxonomy tables/figures, input and output hashes, provenance and SUCCESS.
Fails on upstream inconsistencies, altered reference/inputs, mixed studies or R failure.
On every invocation the prior SUCCESS and validation marker are invalidated first.
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/dada2"))
from validate_outputs import validate_run
from cacao_inventory.config import load_json_yaml
from cacao_inventory.download import validate_file
from cacao_inventory.provenance import file_records
from cacao_inventory.taxonomy_reference import repository_path, validate_taxonomy_config

PRODUCTS = ["taxonomy.tsv", "taxonomy_unfiltered.tsv", "taxonomy_bootstraps.tsv",
            "taxonomy_sensitivity.tsv", "taxonomy_screening.tsv", "assignment_coverage.tsv", "sample_coverage.tsv",
            "software_versions.tsv", "session_info.txt", "warnings.txt", "config_snapshot.yaml",
            "assignment_coverage.pdf", "assignment_coverage.svg", "assignment_coverage.png"]


def write_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/taxonomy.yaml")
    parser.add_argument("--input-dir", type=Path, default=ROOT / "results/dada2/pilot")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/taxonomy/pilot")
    parser.add_argument("--threads", type=int, help="Limit runtime threads to 1..configured threads")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    output = args.output_dir.resolve()
    if ROOT not in output.parents or output == args.input_dir.resolve():
        raise ValueError("Output must be inside the repository and differ from DADA2 inputs")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    status = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).splitlines()
    started = datetime.now(timezone.utc).isoformat()
    output.mkdir(parents=True, exist_ok=True)
    for name in ("SUCCESS", "validation.json"):
        (output / name).unlink(missing_ok=True)
    config = load_json_yaml(args.config)
    validate_taxonomy_config(config, ROOT)
    project = load_json_yaml(ROOT / "config/config.yaml")
    if config["classification"]["random_seed"] != project["project"]["default_random_seed"]:
        raise ValueError("Taxonomy and project seeds disagree")
    threads = args.threads if args.threads is not None else config["classification"]["threads"]
    if not 1 <= threads <= config["classification"]["threads"]:
        raise ValueError("Threads outside the registered range")
    upstream = validate_run(args.input_dir, ROOT)
    reference = config["reference"]
    path = repository_path(ROOT, reference["path"])
    validate_file(path, reference["expected_bytes"], reference["expected_md5"])
    reference_sidecar = repository_path(ROOT, reference["provenance_path"])
    ref_record = json.loads(reference_sidecar.read_text())
    if ref_record["outputs"] != file_records([path], ROOT) or ref_record["parameters"]["reference"] != reference:
        raise ValueError("Reference provenance disagrees with configured reference")
    inputs = [args.config, ROOT / "config/config.yaml", path, reference_sidecar,
              *[args.input_dir / name for name in ("SUCCESS", "validation.json", "asv_sequences.tsv", "asv_counts.tsv", "sample_metadata.tsv")],
              *[ROOT / "scripts/taxonomy" / name for name in ("run_taxonomy.py", "assign_taxonomy.R", "helpers.R")],
              ROOT / "scripts/dada2/helpers.R", ROOT / "scripts/dada2/validate_outputs.py",
              ROOT / "python/cacao_inventory/taxonomy_reference.py", ROOT / "python/cacao_inventory/provenance.py",
              ROOT / "python/cacao_inventory/download.py", ROOT / "python/cacao_inventory/config.py",
              ROOT / "environment/conda-linux-64.lock", ROOT / "workflow/Snakefile", ROOT / "workflow/rules/pilot_taxonomy.smk"]
    before = file_records(inputs, ROOT)
    command = ["Rscript", str(ROOT / "scripts/taxonomy/assign_taxonomy.R"), "--config", str(args.config.resolve()),
               "--input-dir", str(args.input_dir.resolve()), "--output-dir", str(output), "--threads", str(threads)]
    subprocess.run(command, cwd=ROOT, check=True)
    if file_records(inputs, ROOT) != before:
        raise ValueError("Inputs changed during classification; result cannot be marked successful")
    shutil.copyfile(args.config, output / "config_snapshot.yaml")
    write_json(output / "input_checksums.json", before)
    outputs = file_records([output / name for name in PRODUCTS] + [output / "input_checksums.json"], ROOT)
    write_json(output / "provenance.json", {"status": "success", "started_at_utc": started,
        "finished_at_utc": datetime.now(timezone.utc).isoformat(), "git_commit": revision,
        "git_dirty": bool(status), "git_status_at_start": status, "command": command,
        "classification": config["classification"], "rng_kind": ["Mersenne-Twister", "Inversion", "Rejection"],
        "effective_threads": threads, "reference": reference, "upstream_validation": upstream,
        "checksum_algorithm": "SHA-256", "inputs": before, "outputs": outputs})
    (output / "SUCCESS").write_text(revision + "\n", encoding="ascii")
    logging.info("Taxonomy completed for %s ASVs; no exclusions", upstream["asvs"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
