#!/usr/bin/env python3
"""Build the public NCBI/ENA dataset inventory without downloading FASTQ.

Inputs: config/config.yaml and config/datasets.yaml.
Outputs: metadata/studies.tsv, runs.tsv, samples.tsv and exclusion_log.tsv.
Fails on source errors, schema changes, duplicate identifiers or invalid derived
temporal metadata.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_project_configuration  # noqa: E402
from cacao_inventory.inventory import build_inventory  # noqa: E402
from cacao_inventory.io import write_inventory  # noqa: E402
from cacao_inventory.sources import MetadataClient  # noqa: E402
from cacao_inventory.validation import validate_inventory  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "config.yaml")
    parser.add_argument("--datasets", type=Path, default=ROOT / "config" / "datasets.yaml")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "metadata")
    parser.add_argument(
        "--retrieved-at",
        default="",
        help="Optional ISO-8601 UTC timestamp for deterministic rebuilds/tests",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    project_config, candidates = load_project_configuration(args.config, args.datasets)
    client = MetadataClient(project_config["inventory"])
    accessions = [str(candidate["bioproject"]) for candidate in candidates]
    logging.info("Fetching NCBI summaries for %d BioProjects", len(accessions))
    summaries = client.fetch_bioproject_summaries(accessions)
    runs_by_project: dict[str, list[dict[str, str]]] = {}
    for accession in accessions:
        logging.info("Fetching ENA run report for %s", accession)
        runs_by_project[accession] = client.fetch_ena_runs(accession)
    retrieved_at = args.retrieved_at or datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    tables = build_inventory(project_config, candidates, summaries, runs_by_project, retrieved_at)
    summary = validate_inventory(tables, project_config)
    write_inventory(args.output_dir, tables)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
