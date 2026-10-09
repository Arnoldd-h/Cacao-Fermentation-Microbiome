#!/usr/bin/env python3
"""Report immutable local FASTQ validation and public access for missing streams.

Inputs: registered scope/manifest and existing raw files; HEAD requests only.
Outputs: download_status.json and SHA-256 provenance. No downloads or exclusions.
Exit 0 means all streams validate; exit 2 means incomplete/unavailable/altered data.
"""

import argparse
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.analysis_artifacts import write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.download import expand_pilot_manifest, validate_file
from cacao_inventory.fastq import validate_fastq_file
from cacao_inventory.pilot_workflow import load_pilot_manifest
from cacao_inventory.provenance import write_provenance
from cacao_inventory.study_scope import validate_scope


def audit_streams(records, root, opener=urlopen):
    results = []
    valid_paths = []
    for record in records:
        path = root / "data/raw" / record["target_relative"]
        result = dict(record)
        if path.exists():
            try:
                size, checksum = validate_file(path, int(record["expected_bytes"]), record["expected_md5"])
                metrics = validate_fastq_file(path)
                result.update(status="valid", observed_bytes=size, observed_md5=checksum, gzip_valid=True, fastq_valid=True, **metrics)
                valid_paths.append(path)
            except (OSError, EOFError, ValueError) as exc:
                result.update(status="invalid_local_file", error=str(exc))
        else:
            try:
                with opener(Request(record["source_url"], method="HEAD"), timeout=20) as response:
                    kind = response.headers.get("Content-Type", "")
                    result.update(status="html_instead_of_fastq" if kind.lower().startswith("text/html") else "not_downloaded",
                                  http_status=response.status, content_type=kind,
                                  content_length=response.headers.get("Content-Length"), effective_url=response.geturl())
            except Exception as exc:
                result.update(status="remote_request_failed", error=type(exc).__name__ + ": " + str(exc))
        results.append(result)
    valid = sum(row["status"] == "valid" for row in results)
    return {"status": "complete" if valid == len(records) else "blocked", "checked_at_utc": datetime.now(timezone.utc).isoformat(),
            "expected_files": len(records), "valid_files": valid, "pending_files": len(records) - valid,
            "valid_compressed_bytes": sum(row.get("observed_bytes", 0) for row in results),
            "processing_exclusions_added": 0, "streams": results}, valid_paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/full_study.yaml")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    cfg = load_json_yaml(args.config)
    validate_scope(cfg, ROOT)
    manifest = ROOT / cfg["manifest"]
    report, valid_paths = audit_streams(expand_pilot_manifest(load_pilot_manifest(manifest)), ROOT)
    output = ROOT / cfg["plan_dir"] / "download_status.json"
    write_json(output, report)
    write_provenance(output.with_suffix(".provenance.json"), root=ROOT,
                     inputs=[manifest, Path(__file__), ROOT / "python/cacao_inventory/download.py", *valid_paths],
                     outputs=[output], config_paths=[args.config], software_versions={"Python": sys.version.split()[0]}, command=sys.argv)
    logging.info("Validated %d/%d FASTQ files; %d pending", report["valid_files"], report["expected_files"], report["pending_files"])
    return 0 if report["status"] == "complete" else 2


if __name__ == "__main__":
    raise SystemExit(main())
