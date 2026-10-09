#!/usr/bin/env python3
"""Build/audit an exact full-study manifest without downloading sequencing data.

Inputs: scope config, versioned inventory, candidate rules and NCBI/ENA metadata.
Outputs: manifest, source audit, resource plan and hash-based provenance.
Fails on source discrepancies, invalid scope/design/download metadata or low disk.
Source verification is mandatory; pure manifest helpers support offline tests.
"""

import argparse
import logging
import platform
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.download import expand_pilot_manifest, validate_file
from cacao_inventory.io import read_inventory, write_tsv_atomic
from cacao_inventory.pilot_workflow import load_pilot_manifest
from cacao_inventory.provenance import file_records
from cacao_inventory.sources import MetadataClient
from cacao_inventory.study_scope import MANIFEST_COLUMNS, build_study_manifest, build_stream_log, compare_live_runs, validate_scope
from cacao_inventory.validation import validate_inventory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/full_study.yaml")
    parser.add_argument("--verify-sources", action="store_true", required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    cfg = load_json_yaml(args.config)
    validate_scope(cfg, ROOT)
    report = ROOT / cfg["plan_dir"]
    state = begin_stage(report, ROOT)
    inputs = [args.config, ROOT / cfg["project_config"], ROOT / cfg["datasets_config"],
              *[ROOT / cfg["inventory_dir"] / (name + ".tsv") for name in ("studies", "runs", "samples", "exclusion_log")],
              Path(__file__), ROOT / "python/cacao_inventory/study_scope.py", ROOT / "python/cacao_inventory/sources.py"]
    inputs += list((ROOT / "python/cacao_inventory").glob("*.py"))
    before = file_records(inputs, ROOT)
    project = load_json_yaml(ROOT / cfg["project_config"])
    inventory = read_inventory(ROOT / cfg["inventory_dir"])
    validate_inventory(inventory, project)
    rows = build_study_manifest(cfg, inventory["runs"], inventory["samples"])
    stream_log = build_stream_log(cfg, inventory["runs"], inventory["samples"])
    candidates = load_json_yaml(ROOT / cfg["datasets_config"])["candidates"]
    candidate = next(row for row in candidates if row["bioproject"] == cfg["bioproject"])
    if candidate["screening_status"] != "include":
        raise ValueError("Study is not included")
    logging.info("Verifying %s against NCBI and ENA", cfg["bioproject"])
    client = MetadataClient(project["inventory"])
    summary = client.fetch_bioproject_summaries([cfg["bioproject"]])[cfg["bioproject"]]
    audit = compare_live_runs(cfg, candidate, inventory["runs"], client.fetch_ena_runs(cfg["bioproject"]))
    audit["ncbi_bioproject"] = summary["project_acc"]
    audit["ncbi_title"] = summary.get("project_title", "")
    records = expand_pilot_manifest(rows)
    existing = []
    for record in records:
        raw = ROOT / "data/raw" / record["target_relative"]
        if raw.exists():
            validate_file(raw, int(record["expected_bytes"]), record["expected_md5"])
            existing.append(record)
    total = sum(int(row["expected_bytes"]) for row in records)
    required = int(total * cfg["resources"]["minimum_free_space_factor"])
    free = shutil.disk_usage(ROOT).free
    if free < required:
        raise OSError("Insufficient disk for the registered processing reserve")
    plan = {"study_id": cfg["study_id"], "inventory_runs": cfg["expected_runs"], "runs": len(rows), "fastq_files": len(records),
            "excluded_runs_without_paired_files": cfg["expected_runs"] - len(rows),
            "excluded_streams": sum(row["decision"] == "exclude" for row in stream_log),
            "compressed_bytes": total, "existing_valid_files": len(existing),
            "remaining_download_bytes": total - sum(int(row["expected_bytes"]) for row in existing),
            "free_disk_bytes_at_planning": free, "reserved_processing_bytes": required,
            "fermentation_batches": len({row["fermentation_batch"] for row in rows}),
            "batch_time_observations": len({(row["fermentation_batch"], row["fermentation_hours"]) for row in rows}),
            "sampling_strata": sorted({row["sampling_stratum"] for row in rows}), "resources": cfg["resources"],
            "reserve_is_estimate": True, "temporal_inference": False}
    manifest = ROOT / cfg["manifest"]
    write_tsv_atomic(manifest, rows, MANIFEST_COLUMNS)
    load_pilot_manifest(manifest)
    write_json(report / "source_audit.json", audit)
    write_json(report / "resource_plan.json", plan)
    write_tsv_atomic(report / "fastq_stream_log.tsv", stream_log, list(stream_log[0]))
    finish_stage(report, ROOT, state, before, args.config, cfg,
                 [manifest, report / "source_audit.json", report / "resource_plan.json", report / "fastq_stream_log.tsv"],
                 audit, {"Python": platform.python_version()}, sys.argv)
    logging.info("Prepared %d runs, %d FASTQ files, %d compressed bytes; %d valid files reused", len(rows), len(records), total, len(existing))


if __name__ == "__main__":
    main()
