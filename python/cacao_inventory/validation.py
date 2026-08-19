"""Scientific and relational validation for inventory tables."""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from typing import Any

from .config import assign_temporal_stage


class InventoryValidationError(ValueError):
    """Raised when versioned metadata violates a reproducibility invariant."""


def _duplicates(values: list[str]) -> list[str]:
    counts = Counter(values)
    return sorted(value for value, count in counts.items() if value and count > 1)


def validate_inventory(
    tables: dict[str, list[dict[str, str]]], project_config: dict[str, Any]
) -> dict[str, Any]:
    """Validate schemas, keys, temporal derivations and screening invariants."""

    errors: list[str] = []
    studies = tables["studies"]
    runs = tables["runs"]
    samples = tables["samples"]
    exclusions = tables["exclusion_log"]
    stages = project_config["temporal_stages"]
    allowed_decisions = {"true", "false", "pending"}

    duplicate_studies = _duplicates([row.get("study_id", "") for row in studies])
    duplicate_projects = _duplicates([row.get("bioproject", "") for row in studies])
    duplicate_runs = _duplicates([row.get("run_accession", "") for row in runs])
    duplicate_sample_runs = _duplicates([row.get("run_accession", "") for row in samples])
    if duplicate_studies:
        errors.append(f"Duplicate study_id values: {', '.join(duplicate_studies)}")
    if duplicate_projects:
        errors.append(f"Duplicate BioProjects: {', '.join(duplicate_projects)}")
    if duplicate_runs:
        errors.append(f"Duplicate run accessions: {', '.join(duplicate_runs)}")
    if duplicate_sample_runs:
        errors.append(f"Sample rows repeat run accessions: {', '.join(duplicate_sample_runs)}")

    study_ids = {row.get("study_id", "") for row in studies}
    project_ids = {row.get("bioproject", "") for row in studies}
    run_ids = {row.get("run_accession", "") for row in runs}
    for study in studies:
        accession = study.get("bioproject", "")
        if not re.fullmatch(r"PRJ[A-Z]{2}\d+", accession):
            errors.append(f"Invalid BioProject accession: {accession}")
        if study.get("include") not in allowed_decisions:
            errors.append(f"Invalid study include status for {accession}: {study.get('include')}")

    for run in runs:
        accession = run.get("run_accession", "")
        if not re.fullmatch(r"[SED]RR\d+", accession):
            errors.append(f"Invalid run accession: {accession}")
        if run.get("study_id") not in study_ids or run.get("bioproject") not in project_ids:
            errors.append(f"Orphan run: {accession}")
        if run.get("analysis_include") not in allowed_decisions:
            errors.append(f"Invalid run decision: {accession}")

    samples_by_study: dict[str, list[dict[str, str]]] = defaultdict(list)
    for sample in samples:
        run_accession = sample.get("run_accession", "")
        if run_accession not in run_ids:
            errors.append(f"Sample references unknown run: {run_accession}")
        if sample.get("study_id") not in study_ids:
            errors.append(f"Sample references unknown study: {sample.get('study_id')}")
        if sample.get("analysis_include") not in {"true", "pending"}:
            errors.append(f"samples.tsv must contain only true/pending candidates: {run_accession}")
        samples_by_study[sample.get("study_id", "")].append(sample)

        hours_text = sample.get("fermentation_hours", "")
        duration_text = sample.get("fermentation_duration_hours", "")
        relative_text = sample.get("relative_time", "")
        stage = sample.get("fermentation_stage", "")
        if relative_text:
            try:
                hours = float(hours_text)
                duration = float(duration_text)
                relative = float(relative_text)
            except ValueError:
                errors.append(f"Non-numeric temporal metadata: {run_accession}")
                continue
            if duration <= 0 or hours < 0 or hours > duration:
                errors.append(f"Invalid temporal order: {run_accession}")
            expected = hours / duration
            if not math.isclose(relative, expected, rel_tol=1e-7, abs_tol=1e-7):
                errors.append(f"Incorrect relative_time: {run_accession}")
            if not (0 <= relative <= 1):
                errors.append(f"relative_time outside [0,1]: {run_accession}")
            if stage != assign_temporal_stage(relative, stages):
                errors.append(f"Incorrect fermentation_stage: {run_accession}")
        elif stage:
            errors.append(f"Stage without relative_time: {run_accession}")

    for study in studies:
        if study.get("include") != "true":
            continue
        study_samples = samples_by_study.get(study.get("study_id", ""), [])
        included = [row for row in study_samples if row.get("analysis_include") == "true"]
        unique_times = {row.get("fermentation_hours", "") for row in included if row.get("fermentation_hours")}
        if len(unique_times) < 3:
            errors.append(f"Included study lacks at least three time points: {study.get('study_id')}")
        if "16S" not in study.get("marker", ""):
            errors.append(f"Included study is not identified as 16S: {study.get('study_id')}")
        if study.get("raw_data_available") != "true":
            errors.append(f"Included study has no public FASTQ: {study.get('study_id')}")

    logged_run_exclusions = {
        row.get("entity_id", "")
        for row in exclusions
        if row.get("entity_type") == "run" and row.get("decision") == "exclude"
    }
    for run in runs:
        if run.get("analysis_include") == "false" and run.get("run_accession") not in logged_run_exclusions:
            errors.append(f"Excluded run is missing from exclusion_log.tsv: {run.get('run_accession')}")

    if errors:
        preview = "\n- ".join(errors[:30])
        suffix = f"\n... and {len(errors) - 30} more" if len(errors) > 30 else ""
        raise InventoryValidationError(f"Inventory validation failed:\n- {preview}{suffix}")

    decisions = Counter(row.get("include", "") for row in studies)
    run_decisions = Counter(row.get("analysis_include", "") for row in runs)
    return {
        "studies": len(studies),
        "runs": len(runs),
        "samples": len(samples),
        "exclusion_records": len(exclusions),
        "study_decisions": dict(sorted(decisions.items())),
        "run_decisions": dict(sorted(run_decisions.items())),
    }
