#!/usr/bin/env python3
"""Discover cacao-fermentation BioProjects systematically from NCBI and ENA.

This is a broad, reproducible triage. It never promotes a discovered project to
scientific inclusion; priority results still require paper-level screening and a
documented entry in config/datasets.yaml.
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
from cacao_inventory.discovery import build_discovery_rows  # noqa: E402
from cacao_inventory.io import write_tsv_atomic  # noqa: E402
from cacao_inventory.schema import DISCOVERY_COLUMNS  # noqa: E402
from cacao_inventory.sources import MetadataClient  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "config.yaml")
    parser.add_argument("--datasets", type=Path, default=ROOT / "config" / "datasets.yaml")
    parser.add_argument("--output", type=Path, default=ROOT / "metadata" / "discovery_candidates.tsv")
    parser.add_argument("--retrieved-at", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    project_config, candidates = load_project_configuration(args.config, args.datasets)
    client = MetadataClient(project_config["inventory"])
    queries = project_config.get("discovery_queries", [])
    if len(queries) != 1:
        raise ValueError("Bootstrap discovery currently requires exactly one documented query")
    query = queries[0]
    summaries = client.search_bioproject_summaries(str(query["term"]), int(query["retmax"]))
    logging.info("NCBI discovery returned %d BioProjects", len(summaries))
    runs_by_project: dict[str, list[dict[str, str]]] = {}
    for index, accession in enumerate(sorted(summaries), start=1):
        logging.info("Fetching ENA report %d/%d: %s", index, len(summaries), accession)
        runs_by_project[accession] = client.fetch_ena_runs(accession)
    retrieved_at = args.retrieved_at or datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    rows = build_discovery_rows(
        str(query["search_id"]),
        str(query["term"]),
        summaries,
        runs_by_project,
        {str(candidate["bioproject"]) for candidate in candidates},
        retrieved_at,
    )
    write_tsv_atomic(args.output, rows, DISCOVERY_COLUMNS)
    counts: dict[str, int] = {}
    for row in rows:
        key = row["preliminary_screen"]
        counts[key] = counts.get(key, 0) + 1
    print(json.dumps({"discovered": len(rows), "preliminary_screen": counts}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
