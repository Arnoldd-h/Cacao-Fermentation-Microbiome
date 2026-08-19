"""Configuration loading and temporal-stage rules."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


class ConfigurationError(ValueError):
    """Raised when project configuration is missing or inconsistent."""


def load_json_yaml(path: str | Path) -> dict[str, Any]:
    """Load a JSON document stored with a YAML extension.

    JSON is a strict subset of YAML 1.2. Using it here keeps the bootstrap
    inventory dependency-free while remaining consumable by Snakemake/PyYAML.
    """

    source = Path(path)
    try:
        with source.open(encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigurationError(f"Cannot load configuration {source}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigurationError(f"Configuration root must be an object: {source}")
    return data


def validate_temporal_stages(stages: dict[str, Any]) -> None:
    """Validate contiguous early/mid/late intervals spanning zero to one."""

    expected = ("early", "mid", "late")
    intervals: list[tuple[float, float]] = []
    for label in expected:
        value = stages.get(label)
        if not isinstance(value, list) or len(value) != 2:
            raise ConfigurationError(f"temporal_stages.{label} must contain two numbers")
        try:
            lower, upper = float(value[0]), float(value[1])
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(f"temporal_stages.{label} is not numeric") from exc
        if not (math.isfinite(lower) and math.isfinite(upper) and lower <= upper):
            raise ConfigurationError(f"Invalid interval for temporal_stages.{label}")
        intervals.append((lower, upper))

    if not math.isclose(intervals[0][0], 0.0):
        raise ConfigurationError("Temporal stages must start at 0")
    if not math.isclose(intervals[-1][1], 1.0):
        raise ConfigurationError("Temporal stages must end at 1")
    if not math.isclose(intervals[0][1], intervals[1][0]):
        raise ConfigurationError("early and mid boundaries must be contiguous")
    if not math.isclose(intervals[1][1], intervals[2][0]):
        raise ConfigurationError("mid and late boundaries must be contiguous")


def assign_temporal_stage(relative_time: float | None, stages: dict[str, Any]) -> str:
    """Assign early/mid/late using the documented right-closed boundaries."""

    validate_temporal_stages(stages)
    if relative_time is None:
        return ""
    value = float(relative_time)
    if not math.isfinite(value) or value < 0 or value > 1:
        raise ValueError(f"relative_time must be within [0, 1], received {value}")
    early_upper = float(stages["early"][1])
    mid_upper = float(stages["mid"][1])
    if value <= early_upper:
        return "early"
    if value <= mid_upper:
        return "mid"
    return "late"


def load_project_configuration(
    config_path: str | Path, datasets_path: str | Path
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Load and validate project and candidate dataset configuration."""

    project_config = load_json_yaml(config_path)
    datasets_config = load_json_yaml(datasets_path)
    stages = project_config.get("temporal_stages")
    if not isinstance(stages, dict):
        raise ConfigurationError("config.yaml must define temporal_stages")
    validate_temporal_stages(stages)

    candidates = datasets_config.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise ConfigurationError("datasets.yaml must define at least one candidate")

    accessions: set[str] = set()
    study_ids: set[str] = set()
    allowed_statuses = {"include", "exclude", "pending"}
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise ConfigurationError("Each candidate must be an object")
        accession = str(candidate.get("bioproject", ""))
        study_id = str(candidate.get("study_id", ""))
        status = str(candidate.get("screening_status", ""))
        if not accession.startswith("PRJ"):
            raise ConfigurationError(f"Invalid BioProject accession: {accession}")
        if not study_id:
            raise ConfigurationError(f"Missing study_id for {accession}")
        if status not in allowed_statuses:
            raise ConfigurationError(f"Invalid screening_status for {accession}: {status}")
        if accession in accessions:
            raise ConfigurationError(f"Duplicate BioProject accession: {accession}")
        if study_id in study_ids:
            raise ConfigurationError(f"Duplicate study_id: {study_id}")
        accessions.add(accession)
        study_ids.add(study_id)
    return project_config, candidates
