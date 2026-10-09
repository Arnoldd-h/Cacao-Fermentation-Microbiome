"""Exact, inclusive single-study manifests without changing inventory decisions."""

import re
from pathlib import Path

from .download import expand_pilot_manifest
from .pilot_workflow import load_pilot_manifest
from .schema import PILOT_MANIFEST_COLUMNS, SAMPLE_COLUMNS
from .taxonomy_reference import repository_path

MANIFEST_COLUMNS = PILOT_MANIFEST_COLUMNS + [name for name in SAMPLE_COLUMNS if name not in PILOT_MANIFEST_COLUMNS] + ["source_fastq_ftp", "source_fastq_md5", "source_fastq_bytes"]
LIVE_FIELDS = ("sample_alias", "library_strategy", "library_layout", "fastq_ftp", "fastq_bytes", "fastq_md5", "read_count", "base_count")


def validate_scope(config: dict, root: Path) -> None:
    if config["scope"] != "all_included_runs" or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", config["study_id"]):
        raise ValueError("Unsupported scope or unsafe study ID")
    if not re.fullmatch(r"PRJ(?:NA|EB|DB)[0-9]+", config["bioproject"]):
        raise ValueError("Invalid BioProject")
    for field in ("expected_runs", "expected_paired_runs", "expected_batches"):
        if type(config[field]) is not int or config[field] < 1:
            raise ValueError("Expected scope counts must be positive integers")
    if config["expected_paired_runs"] > config["expected_runs"] or config["unpaired_fastq_policy"] != "record_and_exclude_from_paired_processing":
        raise ValueError("Unsupported paired-stream eligibility policy")
    for field in ("manifest", "plan_dir", "qc_dir", "interim_dir", "dada2_dir", "dada2_intermediate_dir", "taxonomy_dir", "project_config", "datasets_config", "inventory_dir", "taxonomy_config", "diversity_config"):
        repository_path(root, config[field])
    for field, parent in (("interim_dir", "data/interim"), ("dada2_intermediate_dir", "results/intermediate"),
                          ("qc_dir", "results/qc"), ("dada2_dir", "results/dada2"), ("taxonomy_dir", "results/taxonomy")):
        if not repository_path(root, config[field]).is_relative_to(root.resolve() / parent) or "/full" not in config[field]:
            raise ValueError("Full-study outputs must use their own scoped directories")
    if config["resources"]["minimum_free_space_factor"] < 2:
        raise ValueError("Insufficient operational disk reserve")
    for field in ("cores", "memory_budget_mb", "fastqc_memory_mb", "analysis_memory_mb"):
        if type(config["resources"][field]) is not int or config["resources"][field] < 1:
            raise ValueError("Invalid execution resource budget")


def select_streams(run: dict, sample_id: str) -> tuple[list[dict], list[dict]]:
    row = {**run, "sample_id": sample_id,
           "fastq_ftp": ";".join(url if url.startswith(("https://", "http://", "ftp://")) else "https://" + url for url in run["fastq_ftp"].split(";"))}
    records = expand_pilot_manifest([row])
    names = {run["run_accession"] + suffix for suffix in ("_1.fastq.gz", "_2.fastq.gz")}
    allowed = names | {run["run_accession"] + ".fastq.gz"}
    if any(record["read_file"] not in allowed or not re.fullmatch(r"[a-fA-F0-9]{32}", record["expected_md5"]) for record in records):
        raise ValueError("Unsupported FASTQ filename or MD5")
    complete = names.issubset({record["read_file"] for record in records})
    paired = [record for record in records if complete and record["read_file"] in names]
    log = [{**record, "decision": "include" if record in paired else "exclude",
            "reason": "paired_stream" if record in paired else "unpaired_stream" if complete else "no_complete_paired_files"} for record in records]
    return sorted(paired, key=lambda record: record["read_file"]), log


def build_stream_log(config: dict, runs: list[dict], samples: list[dict]) -> list[dict]:
    metadata = {row["run_accession"]: row for row in samples if row["study_id"] == config["study_id"] and row["analysis_include"] == "true"}
    logs = []
    for run in runs:
        if run["study_id"] == config["study_id"] and run["analysis_include"] == "true":
            logs.extend(select_streams(run, metadata[run["run_accession"]]["sample_id"])[1])
    return sorted(logs, key=lambda row: (row["run_accession"], row["read_file"]))


def build_study_manifest(config: dict, runs: list[dict], samples: list[dict]) -> list[dict]:
    selected = [row for row in runs if row["study_id"] == config["study_id"] and row["analysis_include"] == "true"]
    metadata = [row for row in samples if row["study_id"] == config["study_id"] and row["analysis_include"] == "true"]
    by_run = {row["run_accession"]: row for row in selected}
    by_sample = {row["run_accession"]: row for row in metadata}
    if len(selected) != config["expected_runs"] or len(by_run) != len(selected) or len(by_sample) != len(metadata) or set(by_run) != set(by_sample):
        raise ValueError("Included runs/samples do not match the exact registered scope")
    if len({row["sample_id"] for row in metadata}) != len(metadata) or len({row["fermentation_batch"] for row in metadata}) != config["expected_batches"]:
        raise ValueError("Duplicate samples or unexpected batch count")
    rows = []
    ordered = sorted(metadata, key=lambda row: (row["fermentation_batch"], float(row["fermentation_hours"]), row["sampling_stratum"], row["run_accession"]))
    for order, sample in enumerate(ordered, 1):
        run = by_run[sample["run_accession"]]
        if sample["bioproject"] != config["bioproject"] or run["bioproject"] != config["bioproject"] or run["library_layout"] != "PAIRED" or run["sequencing_platform"] != "ILLUMINA":
            raise ValueError("Full-study route requires one configured paired Illumina study")
        if not sample["fermentation_batch"] or not sample["sampling_stratum"] or sample["marker"] != "16S rRNA" or not sample["relative_time"]:
            raise ValueError("Missing verified marker, batch, time or stratum")
        records, _ = select_streams(run, sample["sample_id"])
        if not records:
            continue
        row = {**sample, "pilot_order": str(order), "library_layout": run["library_layout"],
               "fastq_ftp": ";".join(record["source_url"] for record in records), "fastq_md5": ";".join(record["expected_md5"] for record in records),
               "fastq_bytes": ";".join(record["expected_bytes"] for record in records), "source_fastq_ftp": run["fastq_ftp"],
               "source_fastq_md5": run["fastq_md5"], "source_fastq_bytes": run["fastq_bytes"],
               "selection_reason": "all inventory-included runs with complete paired streams; no outcome-based subset"}
        row["estimated_bytes_total"] = str(sum(int(r["expected_bytes"]) for r in records))
        rows.append({field: row.get(field, "") for field in MANIFEST_COLUMNS})
    if len(rows) != config["expected_paired_runs"]:
        raise ValueError("Paired processing eligibility changed; review and register the decision")
    if len({row["fermentation_batch"] for row in rows}) != config["expected_batches"]:
        raise ValueError("Paired processing would remove a registered batch")
    for order, row in enumerate(rows, 1):
        row["pilot_order"] = str(order)
    return rows


def compare_live_runs(config: dict, candidate: dict, current: list[dict], live: list[dict]) -> dict:
    live_by_id = {row["run_accession"]: row for row in live}
    if len(live_by_id) != len(live):
        raise ValueError("Duplicate run in ENA response")
    eligible = {row["run_accession"] for row in live if all(re.search(rule["pattern"], row[rule["field"]]) for rule in candidate["candidate_rules"])}
    selected = {row["run_accession"]: row for row in current if row["study_id"] == config["study_id"] and row["analysis_include"] == "true"}
    if eligible != set(selected):
        raise ValueError("ENA eligible run membership changed; review inventory before downloading")
    for accession, original in selected.items():
        if any(original[field] != live_by_id[accession][field] for field in LIVE_FIELDS):
            raise ValueError("ENA download/design metadata changed for " + accession)
    return {"status": "valid", "ena_runs_returned": len(live), "eligible_runs": len(eligible),
            "verified_fields": list(LIVE_FIELDS), "verified_run_accessions": sorted(eligible)}
