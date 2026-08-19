"""Deterministic pilot-dataset comparison from versioned inventory tables."""

from __future__ import annotations

from collections import defaultdict
from typing import Any


TRUE_VALUES = {"true", "pending"}
UNKNOWN_VALUES = {"", "unknown", "not applicable", "not reported"}
STAGE_ORDER = ("early", "mid", "late")


def _is_known(value: str) -> bool:
    normalized = value.strip().lower()
    return normalized not in UNKNOWN_VALUES and not normalized.startswith("unknown")


def _sum_semicolon_integers(value: str) -> int:
    total = 0
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        try:
            total += int(item)
        except ValueError as exc:
            raise ValueError(f"Invalid FASTQ byte count: {item}") from exc
    return total


def _format_hours(values: set[float]) -> str:
    return ";".join(f"{value:g}" for value in sorted(values))


def _ranking_key(row: dict[str, str], rules: dict[str, Any]) -> tuple[Any, ...]:
    quality_order = [str(value) for value in rules["metadata_quality_order"]]
    quality_rank = {label: index for index, label in enumerate(quality_order)}
    key_parts: list[Any] = []
    for field in rules["ranking_priority"]:
        if field == "metadata_quality":
            key_parts.append(quality_rank.get(row[field], len(quality_order)))
        elif field in {"primers_known", "paired_end"}:
            key_parts.append(0 if row[field] == "true" else 1)
        elif field in {"unique_timepoints", "fermentation_batches"}:
            key_parts.append(-int(row[field]))
        elif field == "selected_fastq_bytes":
            key_parts.append(int(row[field]))
        else:
            raise ValueError(f"Unsupported pilot ranking field: {field}")
    key_parts.append(row["study_id"])
    return tuple(key_parts)


def build_pilot_selection(
    studies: list[dict[str, str]],
    runs: list[dict[str, str]],
    samples: list[dict[str, str]],
    rules: dict[str, Any],
) -> list[dict[str, str]]:
    """Compare all studies and identify the highest-ranked eligible primary study."""

    minimum_timepoints = int(rules["minimum_timepoints"])
    required_stages = {str(value) for value in rules["required_temporal_stages"]}
    if minimum_timepoints < 3:
        raise ValueError("pilot_selection.minimum_timepoints must be at least 3")
    if not required_stages:
        raise ValueError("pilot_selection.required_temporal_stages cannot be empty")

    runs_by_study: dict[str, list[dict[str, str]]] = defaultdict(list)
    samples_by_study: dict[str, list[dict[str, str]]] = defaultdict(list)
    for run in runs:
        if run["analysis_include"] in TRUE_VALUES:
            runs_by_study[run["study_id"]].append(run)
    for sample in samples:
        if sample["analysis_include"] in TRUE_VALUES:
            samples_by_study[sample["study_id"]].append(sample)

    rows: list[dict[str, str]] = []
    for study in studies:
        study_id = study["study_id"]
        study_runs = runs_by_study[study_id]
        study_samples = samples_by_study[study_id]
        times = {
            float(sample["fermentation_hours"])
            for sample in study_samples
            if sample["fermentation_hours"]
        }
        stages = {
            sample["fermentation_stage"]
            for sample in study_samples
            if sample["fermentation_stage"]
        }
        batches = {
            sample["fermentation_batch"]
            for sample in study_samples
            if sample["fermentation_batch"]
        }
        selected_bytes = sum(_sum_semicolon_integers(run["fastq_bytes"]) for run in study_runs)
        publication_linked = bool(study["doi"] or study["pmid"])
        primers_known = _is_known(study["forward_primer"]) and _is_known(
            study["reverse_primer"]
        )

        failures: list[str] = []
        if study["include"] != "true":
            failures.append(f"study decision is {study['include'] or 'unknown'}")
        if study["raw_data_available"] != "true":
            failures.append("raw FASTQ is not verified as public")
        if "16s" not in study["marker"].lower():
            failures.append("bacterial 16S marker is not verified")
        if len(times) < minimum_timepoints:
            failures.append(f"fewer than {minimum_timepoints} timepoints")
        missing_stages = sorted(required_stages - stages)
        if missing_stages:
            failures.append(f"missing temporal stages: {','.join(missing_stages)}")

        eligible = not failures
        rows.append(
            {
                "selection_rank": "",
                "selected_as_pilot": "false",
                "eligible_primary": str(eligible).lower(),
                "study_id": study_id,
                "bioproject": study["bioproject"],
                "study_include": study["include"],
                "metadata_quality": study["metadata_quality"],
                "publication_linked": str(publication_linked).lower(),
                "primers_known": str(primers_known).lower(),
                "paired_end": study["paired_end"],
                "raw_data_available": study["raw_data_available"],
                "candidate_runs": str(len(study_runs)),
                "fermentation_batches": str(len(batches)),
                "unique_timepoints": str(len(times)),
                "timepoints_hours": _format_hours(times),
                "temporal_stages": ";".join(stage for stage in STAGE_ORDER if stage in stages),
                "selected_fastq_bytes": str(selected_bytes),
                "selected_fastq_mib": f"{selected_bytes / 1024**2:.2f}",
                "eligibility_reason": "eligible" if eligible else "; ".join(failures),
            }
        )

    eligible_rows = [row for row in rows if row["eligible_primary"] == "true"]
    eligible_rows.sort(key=lambda row: _ranking_key(row, rules))
    for rank, row in enumerate(eligible_rows, start=1):
        row["selection_rank"] = str(rank)
        row["selected_as_pilot"] = str(rank == 1).lower()

    rows.sort(
        key=lambda row: (
            row["eligible_primary"] != "true",
            int(row["selection_rank"]) if row["selection_rank"] else 10**9,
            row["study_id"],
        )
    )
    return rows
