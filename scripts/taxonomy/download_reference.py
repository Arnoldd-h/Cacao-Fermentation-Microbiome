#!/usr/bin/env python3
"""Download the pinned taxonomy reference and record validated provenance.

Input: JSON-compatible taxonomy configuration with URL, bytes and provider MD5.
Outputs: ignored gzip FASTA under references/databases and a small SHA-256 sidecar.
Fails on invalid configuration, insufficient disk, checksum or FASTA inconsistency.
Existing files are validated and never silently replaced. Partial downloads resume.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from functools import partial
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.config import load_json_yaml
from cacao_inventory.download import download_file
from cacao_inventory.provenance import write_provenance
from cacao_inventory.taxonomy_reference import inspect_training_fasta, repository_path, validate_taxonomy_config


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/taxonomy.yaml")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    config = load_json_yaml(args.config)
    validate_taxonomy_config(config, ROOT)
    reference, settings = config["reference"], config["download"]
    path = repository_path(ROOT, reference["path"])
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() and shutil.disk_usage(path.parent).free < reference["expected_bytes"] * settings["minimum_free_space_factor"]:
        raise ValueError("Insufficient disk for reference download and processing")
    for attempt in range(int(settings["retries"])):
        try:
            status = download_file(reference["url"], path, reference["expected_bytes"], reference["expected_md5"],
                                   settings["chunk_size_bytes"], opener=partial(urlopen, timeout=settings["timeout_seconds"]))
            break
        except (URLError, TimeoutError, ConnectionError, OSError):
            if attempt + 1 == settings["retries"]:
                raise
            logging.warning("Network failure; retrying reference download (%s)", attempt + 1)
    summary = inspect_training_fasta(path, len(config["classification"]["tax_levels"]))
    write_provenance(repository_path(ROOT, reference["provenance_path"]), root=ROOT,
                     inputs=[Path(__file__), ROOT / "python/cacao_inventory/taxonomy_reference.py",
                             ROOT / "python/cacao_inventory/download.py"], outputs=[path], config_paths=[args.config],
                     software_versions={"Python": sys.version.split()[0]}, command=sys.argv,
                     parameters={"reference": reference, "fasta_inspection": summary})
    logging.info("Reference %s: %s; %s", reference["version"], status, summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
