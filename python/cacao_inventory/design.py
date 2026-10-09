"""Explicit nesting of sequencing runs within longitudinal fermentation units."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from typing import Any

UNIT_COLUMNS = [
    "study_id", "bioproject", "sample_id", "run_accession", "fermentation_batch",
    "fermentation_hours", "observation_id", "subsamples_in_observation",
    "sampling_stratum", "replicate", "relative_time", "time_source",
]
DESIGN_COLUMNS = [
    "study_id", "bioproject", "country", "region_16s", "runs",
    "independent_batches", "batch_time_observations", "unique_timepoints",
    "observations_with_subsamples", "minimum_subsamples", "maximum_subsamples",
]


def canonical_hours(value: str) -> str:
    """Normalize numerical hours so 24 and 24.0 identify the same observation."""
    try:
        number = Decimal(value)
    except (InvalidOperation, TypeError) as exc:
        raise ValueError(f"Invalid fermentation_hours: {value!r}") from exc
    if not number.is_finite() or number < 0:
        raise ValueError(f"Invalid fermentation_hours: {value!r}")
    return format(number.normalize(), "f")


def build_analysis_units(
    samples: list[dict[str, str]], studies: list[dict[str, str]], config: dict[str, Any]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Return run-level nesting and study summaries without pooling abundances."""
    settings = config["analysis_design"]
    if settings["independent_unit"] != "fermentation_batch":
        raise ValueError("Only the preregistered fermentation_batch unit is implemented")
    fields = settings["observation_fields"]
    if not isinstance(fields, list) or len(fields) != 3 or set(fields) != {"study_id", "fermentation_batch", "fermentation_hours"}:
        raise ValueError("Observation fields must preserve study, batch and time")
    studies_by_id = {row["study_id"]: row for row in studies}
    if len(studies_by_id) != len(studies):
        raise ValueError("Duplicate study_id")
    rows: list[dict[str, str]] = []
    seen_runs: set[str] = set()
    for sample in samples:
        if sample["analysis_include"] != "true":
            continue
        study = studies_by_id.get(sample["study_id"])
        if not study or study["include"] != "true":
            raise ValueError("Included sample does not belong to an included study")
        if sample["bioproject"] != study["bioproject"]:
            raise ValueError("Sample BioProject disagrees with its study")
        if not sample["run_accession"] or sample["run_accession"] in seen_runs:
            raise ValueError("Missing or duplicate run_accession")
        seen_runs.add(sample["run_accession"])
        row = {key: sample.get(key, "") for key in UNIT_COLUMNS}
        row["fermentation_hours"] = canonical_hours(sample["fermentation_hours"])
        if any(not row[field] for field in fields):
            raise ValueError(f"Incomplete longitudinal unit: {sample['run_accession']}")
        identity = json.dumps([row[field] for field in fields], ensure_ascii=True)
        row["observation_id"] = "OBS_" + hashlib.sha256(identity.encode()).hexdigest()[:16]
        rows.append(row)
    counts = Counter(row["observation_id"] for row in rows)
    by_study: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        row["subsamples_in_observation"] = str(counts[row["observation_id"]])
        by_study[row["study_id"]].append(row)
    summaries: list[dict[str, str]] = []
    for study_id, members in sorted(by_study.items()):
        study = studies_by_id[study_id]
        observation_counts = Counter(row["observation_id"] for row in members)
        summaries.append({
            "study_id": study_id, "bioproject": study["bioproject"],
            "country": study["country"], "region_16s": study["region_16s"],
            "runs": str(len(members)),
            "independent_batches": str(len({row['fermentation_batch'] for row in members})),
            "batch_time_observations": str(len(observation_counts)),
            "unique_timepoints": str(len({row['fermentation_hours'] for row in members})),
            "observations_with_subsamples": str(sum(count > 1 for count in observation_counts.values())),
            "minimum_subsamples": str(min(observation_counts.values())),
            "maximum_subsamples": str(max(observation_counts.values())),
        })
    return sorted(rows, key=lambda row: (row["study_id"], row["run_accession"])), summaries
