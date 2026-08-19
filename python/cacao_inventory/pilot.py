"""Build a small, balanced FASTQ manifest for pilot workflow validation."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from typing import Any


def _parts(value: str) -> list[str]:
    return [part.strip() for part in value.split(";") if part.strip()]


def _https_fastq_paths(value: str) -> str:
    paths = []
    for path in _parts(value):
        if path.startswith(("https://", "http://", "ftp://")):
            paths.append(path)
        else:
            paths.append(f"https://{path}")
    return ";".join(paths)


def _validate_download_fields(run: dict[str, str], required_files: int) -> int | None:
    paths = _parts(run["fastq_ftp"])
    checksums = _parts(run["fastq_md5"])
    byte_values = _parts(run["fastq_bytes"])
    if not (len(paths) == len(checksums) == len(byte_values) == required_files):
        return None
    if any(re.fullmatch(r"[0-9a-fA-F]{32}", checksum) is None for checksum in checksums):
        return None
    try:
        parsed_bytes = [int(value) for value in byte_values]
    except ValueError:
        return None
    if any(value <= 0 for value in parsed_bytes):
        return None
    return sum(parsed_bytes)


def build_pilot_manifest(
    pilot_study_id: str,
    runs: list[dict[str, str]],
    samples: list[dict[str, str]],
    rules: dict[str, Any],
) -> list[dict[str, str]]:
    """Select independent batches across temporal stages with complete FASTQ metadata."""

    samples_per_stage = int(rules["samples_per_stage"])
    required_stages = [str(value) for value in rules["required_stages"]]
    required_layout = str(rules["required_layout"])
    required_files = int(rules["required_fastq_files"])
    if samples_per_stage < 1 or required_files < 1:
        raise ValueError("Pilot manifest counts must be positive")
    if len(required_stages) != len(set(required_stages)):
        raise ValueError("pilot_vertical_slice.required_stages must be unique")

    run_by_accession = {
        run["run_accession"]: run
        for run in runs
        if run["study_id"] == pilot_study_id and run["analysis_include"] == "true"
    }
    candidates_by_stage: dict[str, list[tuple[dict[str, str], dict[str, str], int]]] = (
        defaultdict(list)
    )
    for sample in samples:
        if sample["study_id"] != pilot_study_id or sample["analysis_include"] != "true":
            continue
        run = run_by_accession.get(sample["run_accession"])
        if run is None or run["library_layout"] != required_layout:
            continue
        total_bytes = _validate_download_fields(run, required_files)
        if total_bytes is None:
            continue
        stage = sample["fermentation_stage"]
        if stage in required_stages:
            candidates_by_stage[stage].append((sample, run, total_bytes))

    batch_usage: Counter[str] = Counter()
    stratum_usage: Counter[str] = Counter()
    selected: list[tuple[dict[str, str], dict[str, str], int]] = []
    for stage in required_stages:
        stage_candidates = candidates_by_stage[stage]
        available_batches = {
            sample["fermentation_batch"] for sample, _, _ in stage_candidates if sample["fermentation_batch"]
        }
        if len(available_batches) < samples_per_stage:
            raise ValueError(
                f"Stage {stage} has {len(available_batches)} eligible batches; "
                f"{samples_per_stage} required"
            )
        stage_batches: set[str] = set()
        for _ in range(samples_per_stage):
            remaining = [
                candidate
                for candidate in stage_candidates
                if candidate[0]["fermentation_batch"] not in stage_batches
            ]
            remaining.sort(
                key=lambda candidate: (
                    batch_usage[candidate[0]["fermentation_batch"]],
                    stratum_usage[candidate[0]["sampling_stratum"]],
                    float(candidate[0]["fermentation_hours"]),
                    candidate[0]["fermentation_batch"],
                    candidate[0]["sampling_stratum"],
                    candidate[0]["sample_id"],
                    candidate[1]["run_accession"],
                )
            )
            chosen = remaining[0]
            selected.append(chosen)
            batch = chosen[0]["fermentation_batch"]
            stratum = chosen[0]["sampling_stratum"]
            stage_batches.add(batch)
            batch_usage[batch] += 1
            if stratum:
                stratum_usage[stratum] += 1

    rows: list[dict[str, str]] = []
    for order, (sample, run, total_bytes) in enumerate(selected, start=1):
        rows.append(
            {
                "pilot_order": str(order),
                "study_id": pilot_study_id,
                "bioproject": sample["bioproject"],
                "sample_id": sample["sample_id"],
                "run_accession": run["run_accession"],
                "fermentation_stage": sample["fermentation_stage"],
                "fermentation_hours": sample["fermentation_hours"],
                "relative_time": sample["relative_time"],
                "fermentation_batch": sample["fermentation_batch"],
                "sampling_stratum": sample["sampling_stratum"],
                "library_layout": run["library_layout"],
                "sequencing_platform": run["sequencing_platform"],
                "sequencing_instrument": run["sequencing_instrument"],
                "fastq_ftp": _https_fastq_paths(run["fastq_ftp"]),
                "fastq_md5": run["fastq_md5"],
                "fastq_bytes": run["fastq_bytes"],
                "estimated_bytes_total": str(total_bytes),
                "selection_reason": (
                    f"balanced {sample['fermentation_stage']} stage; independent fermentation "
                    "batch; paired FASTQ paths, bytes and MD5 present in ENA metadata"
                ),
            }
        )
    return rows
