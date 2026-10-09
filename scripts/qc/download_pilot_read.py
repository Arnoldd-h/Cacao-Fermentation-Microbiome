#!/usr/bin/env python3
"""Download one manifest FASTQ, checking bytes/MD5 before atomic promotion.

Inputs: pilot manifest, config, study/run/direction. Output: its exact data/raw
path. Fails on ambiguous metadata, invalid paths, insufficient disk, or checksum
mismatch. Existing valid raw files are reused; invalid raw files are never replaced.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_json_yaml
from cacao_inventory.download import download_file, expand_pilot_manifest
from cacao_inventory.pilot_workflow import fastq_path, load_pilot_manifest, manifest_row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--run-accession", required=True)
    parser.add_argument("--direction", choices=("1", "2"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    row = manifest_row(load_pilot_manifest(args.manifest), args.run_accession, args.study_id)
    expected = ROOT / fastq_path(row, args.direction)
    if args.output.resolve() != expected.resolve():
        raise ValueError("Output path disagrees with the immutable manifest raw path")
    record = next(record for record in expand_pilot_manifest([row]) if record["read_file"] == expected.name)
    settings = load_json_yaml(args.config)["pilot_download"]
    chunk_size = int(settings["chunk_size_bytes"])
    factor = float(settings["minimum_free_space_factor"])
    if chunk_size <= 0 or factor < 1:
        raise ValueError("Invalid pilot_download settings")
    expected.parent.mkdir(parents=True, exist_ok=True)
    if not expected.exists() and shutil.disk_usage(expected.parent).free < int(record["expected_bytes"]) * factor:
        raise OSError("Insufficient disk space for pilot FASTQ")
    status = download_file(record["source_url"], expected, int(record["expected_bytes"]), record["expected_md5"], chunk_size)
    logging.info("%s: %s", expected.name, status)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
