"""Transform authoritative source rows into auditable project inventory tables."""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any

from .config import assign_temporal_stage
from .schema import EXCLUSION_COLUMNS, RUN_COLUMNS, SAMPLE_COLUMNS, STUDY_COLUMNS


def _blank(columns: list[str]) -> dict[str, str]:
    return {column: "" for column in columns}


def _regex_context(row: dict[str, str], rules: list[dict[str, str]]) -> dict[str, str] | None:
    context: dict[str, str] = {}
    if not rules:
        return None
    for rule in rules:
        field = rule.get("field", "")
        pattern = rule.get("pattern", "")
        match = re.search(pattern, row.get(field, ""), flags=re.IGNORECASE)
        if match is None:
            return None
        context.update({key: value for key, value in match.groupdict().items() if value is not None})
    return context


def _parse_time(
    row: dict[str, str], candidate: dict[str, Any], base_context: dict[str, str]
) -> tuple[float | None, dict[str, str], str]:
    parser = candidate.get("time_parser")
    if not isinstance(parser, dict):
        return None, dict(base_context), ""
    field = str(parser.get("field", ""))
    match = re.search(str(parser.get("pattern", "")), row.get(field, ""), flags=re.IGNORECASE)
    if match is None:
        return None, dict(base_context), ""
    context = dict(base_context)
    context.update({key: value for key, value in match.groupdict().items() if value is not None})
    try:
        zero_group = str(parser.get("zero_group", ""))
        if zero_group and context.get(zero_group):
            hours = 0.0
        elif parser.get("unit") == "days":
            hours = float(context["days"]) * 24.0
        else:
            hours = float(context["hours"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Time parser for {candidate['bioproject']} did not produce numeric time") from exc
    return hours, context, f"ENA {field} parsed with configured regex"


def _format_number(value: float | int | None) -> str:
    if value is None:
        return ""
    numeric = float(value)
    if numeric.is_integer():
        return str(int(numeric))
    return f"{numeric:.8g}"


def _sum_semicolon_integers(values: list[str]) -> int:
    total = 0
    for value in values:
        for part in value.split(";"):
            stripped = part.strip()
            if stripped:
                try:
                    total += int(stripped)
                except ValueError:
                    continue
    return total


def _split_country(value: str, fallback_country: str, fallback_location: str) -> tuple[str, str]:
    if ":" in value:
        country, location = value.split(":", 1)
        return country.strip(), location.strip()
    if value.strip():
        return value.strip(), fallback_location
    return fallback_country, fallback_location


def _mapped_value(candidate: dict[str, Any], key: str, context: dict[str, str]) -> str:
    raw = context.get(key, "")
    mapping = candidate.get("value_maps", {}).get(key, {})
    return str(mapping.get(raw, raw))


def build_samples(
    candidate: dict[str, Any], rows: list[dict[str, str]], stages: dict[str, Any]
) -> tuple[list[dict[str, str]], set[str]]:
    """Build candidate 16S sample rows using dataset-specific configured rules."""

    selected: list[tuple[dict[str, str], dict[str, str], float | None, str]] = []
    for row in rows:
        context = _regex_context(row, candidate.get("candidate_rules", []))
        if context is None:
            continue
        hours, parsed_context, time_source = _parse_time(row, candidate, context)
        selected.append((row, parsed_context, hours, time_source))

    group_durations: dict[tuple[str, ...], float] = {}
    if candidate.get("duration_mode") == "group_max":
        group_fields = [str(value) for value in candidate.get("duration_group_fields", [])]
        grouped: dict[tuple[str, ...], list[float]] = defaultdict(list)
        for _, context, hours, _ in selected:
            if hours is not None:
                grouped[tuple(context.get(field, "") for field in group_fields)].append(hours)
        group_durations = {key: max(values) for key, values in grouped.items() if values}

    output: list[dict[str, str]] = []
    selected_accessions: set[str] = set()
    status = str(candidate["screening_status"])
    decision = "true" if status == "include" else "pending" if status == "pending" else "false"
    group_fields = [str(value) for value in candidate.get("duration_group_fields", [])]

    for row, context, hours, time_source in selected:
        duration: float | None = None
        if candidate.get("duration_mode") == "fixed":
            duration = float(candidate["duration_hours"])
        elif candidate.get("duration_mode") == "group_max":
            duration = group_durations.get(tuple(context.get(field, "") for field in group_fields))

        relative_time: float | None = None
        if hours is not None and duration is not None:
            if duration <= 0 or hours < 0 or hours > duration:
                raise ValueError(
                    f"Invalid time/duration for {row.get('run_accession')}: {hours}/{duration}"
                )
            relative_time = hours / duration

        raw_country = row.get("country", "")
        country, raw_location = _split_country(
            raw_country, str(candidate.get("country", "")), str(candidate.get("region", ""))
        )
        location = _mapped_value(candidate, "location_code", context) or raw_location
        season = _mapped_value(candidate, "season_code", context) or str(candidate.get("season", ""))
        group_values = [context.get(field, "") for field in group_fields]
        fermentation_batch = str(candidate["study_id"])
        if any(group_values):
            fermentation_batch += "::" + "::".join(group_values)

        sample_variety = _mapped_value(candidate, "variety", context) or str(
            candidate.get("cacao_variety", "")
        )
        sample = _blank(SAMPLE_COLUMNS)
        sample.update(
            {
                "study_id": str(candidate["study_id"]),
                "sample_id": row.get("sample_alias") or row.get("sample_accession", ""),
                "run_accession": row.get("run_accession", ""),
                "bioproject": str(candidate["bioproject"]),
                "biosample": row.get("sample_accession", ""),
                "sample_alias": row.get("sample_alias", ""),
                "country": country,
                "location": location,
                "fermentation_batch": fermentation_batch,
                "fermentation_hours": _format_number(hours),
                "fermentation_duration_hours": _format_number(duration),
                "relative_time": _format_number(relative_time),
                "fermentation_stage": assign_temporal_stage(relative_time, stages),
                "season": season,
                "cacao_variety": sample_variety,
                "replicate": context.get("replicate", ""),
                "sampling_stratum": context.get("sampling_stratum", ""),
                "sequencing_platform": row.get("instrument_platform", ""),
                "sequencing_instrument": row.get("instrument_model", ""),
                "marker": "16S rRNA",
                "region_16s": str(candidate.get("region_16s", "")),
                "time_source": time_source,
                "analysis_include": decision,
                "exclusion_reason": "" if decision == "true" else str(candidate.get("screening_reason", "")),
                "notes": "; ".join(
                    value
                    for value in [row.get("sample_title", ""), row.get("sample_description", "")]
                    if value
                ),
            }
        )
        output.append(sample)
        selected_accessions.add(sample["run_accession"])

    output.sort(
        key=lambda item: (
            item["study_id"],
            item["fermentation_batch"],
            float(item["fermentation_hours"]) if item["fermentation_hours"] else float("inf"),
            item["run_accession"],
        )
    )
    return output, selected_accessions


def _non_candidate_reason(row: dict[str, str], candidate: dict[str, Any]) -> tuple[str, str]:
    text = " ".join(
        row.get(field, "")
        for field in ("sample_alias", "sample_title", "sample_description", "experiment_title")
    ).lower()
    if row.get("library_strategy", "").upper() == "WGS":
        return "whole_metagenome_shotgun", "WGS is outside Phase I bacterial 16S amplicon scope"
    biological_control = bool(re.search(r"fermentation\s+\d+.*negative control", text))
    if "mock" in text or ("negative control" in text and not biological_control):
        return "control", "Mock or negative-control run is not a biological fermentation sample"
    if "phyllosphere" in text:
        return "non_fermentation_matrix", "Phyllosphere run is outside the fermentation-mass population"
    non_bacterial_signal = any(value in text for value in ("its", "fung", "yeast"))
    bacterial_signal = "16s" in text or "bacter" in text
    if non_bacterial_signal and not bacterial_signal:
        return "non_bacterial_marker", "ITS/fungal run is outside bacterial 16S scope"
    if candidate.get("screening_status") == "exclude":
        return "study_excluded", str(candidate.get("screening_reason", "Study excluded"))
    return "not_verified_16s_subset", "Run does not match the configured, verified bacterial 16S subset"


def build_runs(
    candidate: dict[str, Any], rows: list[dict[str, str]], selected_accessions: set[str]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build all-run inventory and explicit run-level exclusions."""

    run_rows: list[dict[str, str]] = []
    exclusions: list[dict[str, str]] = []
    status = str(candidate["screening_status"])
    selected_decision = "true" if status == "include" else "pending" if status == "pending" else "false"
    source = f"https://www.ebi.ac.uk/ena/browser/view/{candidate['bioproject']}"

    for raw in rows:
        selected = raw.get("run_accession", "") in selected_accessions
        if selected:
            marker_classification = "bacterial_16s_candidate"
            decision = selected_decision
            reason = "" if decision == "true" else str(candidate.get("screening_reason", ""))
        else:
            marker_classification, reason = _non_candidate_reason(raw, candidate)
            decision = "false"

        run = _blank(RUN_COLUMNS)
        run.update(
            {
                "study_id": str(candidate["study_id"]),
                "bioproject": str(candidate["bioproject"]),
                "study_accession": raw.get("study_accession", ""),
                "secondary_study_accession": raw.get("secondary_study_accession", ""),
                "run_accession": raw.get("run_accession", ""),
                "experiment_accession": raw.get("experiment_accession", ""),
                "biosample": raw.get("sample_accession", ""),
                "secondary_sample_accession": raw.get("secondary_sample_accession", ""),
                "sample_alias": raw.get("sample_alias", ""),
                "sample_title": raw.get("sample_title", ""),
                "sample_description": raw.get("sample_description", ""),
                "experiment_title": raw.get("experiment_title", ""),
                "scientific_name": raw.get("scientific_name", ""),
                "tax_id": raw.get("tax_id", ""),
                "library_strategy": raw.get("library_strategy", ""),
                "library_source": raw.get("library_source", ""),
                "library_selection": raw.get("library_selection", ""),
                "library_layout": raw.get("library_layout", ""),
                "sequencing_platform": raw.get("instrument_platform", ""),
                "sequencing_instrument": raw.get("instrument_model", ""),
                "base_count": raw.get("base_count", ""),
                "read_count": raw.get("read_count", ""),
                "fastq_bytes": raw.get("fastq_bytes", ""),
                "fastq_ftp": raw.get("fastq_ftp", ""),
                "fastq_md5": raw.get("fastq_md5", ""),
                "country": raw.get("country", ""),
                "collection_date": raw.get("collection_date", ""),
                "cultivar": raw.get("cultivar", ""),
                "isolation_source": raw.get("isolation_source", ""),
                "first_public": raw.get("first_public", ""),
                "last_updated": raw.get("last_updated", ""),
                "marker_classification": marker_classification,
                "analysis_include": decision,
                "exclusion_reason": reason,
                "source": source,
            }
        )
        run_rows.append(run)

        if decision == "false":
            exclusion = _blank(EXCLUSION_COLUMNS)
            exclusion.update(
                {
                    "entity_type": "run",
                    "entity_id": run["run_accession"],
                    "study_id": run["study_id"],
                    "bioproject": run["bioproject"],
                    "reason": reason,
                    "metric": marker_classification,
                    "threshold": "Phase I bacterial 16S fermentation time-series",
                    "decision": "exclude",
                    "source": source,
                }
            )
            exclusions.append(exclusion)
    run_rows.sort(key=lambda item: item["run_accession"])
    exclusions.sort(key=lambda item: item["entity_id"])
    return run_rows, exclusions


def build_study(
    candidate: dict[str, Any],
    summary: dict[str, Any],
    runs: list[dict[str, str]],
    samples: list[dict[str, str]],
    retrieved_at_utc: str,
) -> dict[str, str]:
    """Build one study-level record from deposited and curated evidence."""

    platforms = sorted({row.get("instrument_platform", "") for row in runs if row.get("instrument_platform")})
    instruments = sorted({row.get("instrument_model", "") for row in runs if row.get("instrument_model")})
    layouts = sorted({row.get("library_layout", "") for row in runs if row.get("library_layout")})
    if layouts == ["PAIRED"]:
        paired_end = "true"
    elif layouts == ["SINGLE"]:
        paired_end = "false"
    else:
        paired_end = "mixed" if layouts else ""
    status = str(candidate["screening_status"])
    include = "true" if status == "include" else "false" if status == "exclude" else "pending"

    study = _blank(STUDY_COLUMNS)
    study.update(
        {
            "study_id": str(candidate["study_id"]),
            "bioproject": str(candidate["bioproject"]),
            "publication_title": str(candidate.get("publication_title", "")),
            "doi": str(candidate.get("doi", "")),
            "pmid": str(candidate.get("pmid", "")),
            "year": str(candidate.get("year", "")),
            "project_title": str(summary.get("project_title", "")),
            "project_description": str(summary.get("project_description", "")),
            "country": str(candidate.get("country", "")),
            "region": str(candidate.get("region", "")),
            "institution": str(candidate.get("institution", "")),
            "cacao_species": str(candidate.get("cacao_species", "")),
            "cacao_variety": str(candidate.get("cacao_variety", "")),
            "fermentation_type": str(candidate.get("fermentation_type", "")),
            "fermentation_container": str(candidate.get("fermentation_container", "")),
            "fermentation_duration_hours": str(candidate.get("fermentation_duration_hours", "")),
            "season": str(candidate.get("season", "")),
            "amplicon_type": str(candidate.get("amplicon_type", "")),
            "marker": str(candidate.get("marker", "")),
            "region_16s": str(candidate.get("region_16s", "")),
            "forward_primer": str(candidate.get("forward_primer", "")),
            "reverse_primer": str(candidate.get("reverse_primer", "")),
            "sequencing_platform": "; ".join(platforms),
            "sequencing_instrument": "; ".join(instruments),
            "paired_end": paired_end,
            "total_runs": str(len(runs)),
            "number_samples": str(len(samples)),
            "estimated_fastq_bytes": str(_sum_semicolon_integers([row.get("fastq_bytes", "") for row in runs])),
            "raw_data_available": "true" if any(row.get("fastq_ftp") for row in runs) else "false",
            "metadata_quality": str(candidate.get("metadata_quality", "")),
            "include": include,
            "exclusion_reason": "" if include == "true" else str(candidate.get("screening_reason", "")),
            "source": " | ".join(str(value) for value in candidate.get("sources", [])),
            "retrieved_at_utc": retrieved_at_utc,
        }
    )
    return study


def build_inventory(
    project_config: dict[str, Any],
    candidates: list[dict[str, Any]],
    summaries: dict[str, dict[str, Any]],
    runs_by_project: dict[str, list[dict[str, str]]],
    retrieved_at_utc: str,
) -> dict[str, list[dict[str, str]]]:
    """Build the four versioned inventory tables."""

    studies: list[dict[str, str]] = []
    all_runs: list[dict[str, str]] = []
    all_samples: list[dict[str, str]] = []
    all_exclusions: list[dict[str, str]] = []
    stages = project_config["temporal_stages"]

    for candidate in candidates:
        accession = str(candidate["bioproject"])
        raw_runs = runs_by_project.get(accession, [])
        samples, selected_accessions = build_samples(candidate, raw_runs, stages)
        run_rows, run_exclusions = build_runs(candidate, raw_runs, selected_accessions)
        study = build_study(candidate, summaries[accession], raw_runs, samples, retrieved_at_utc)
        studies.append(study)
        all_runs.extend(run_rows)
        all_samples.extend(samples)
        all_exclusions.extend(run_exclusions)

        if candidate["screening_status"] in {"pending", "exclude"}:
            exclusion = _blank(EXCLUSION_COLUMNS)
            exclusion.update(
                {
                    "entity_type": "study",
                    "entity_id": accession,
                    "study_id": str(candidate["study_id"]),
                    "bioproject": accession,
                    "reason": str(candidate.get("screening_reason", "")),
                    "metric": "screening_status",
                    "threshold": "protocol/inclusion_criteria.md",
                    "decision": str(candidate["screening_status"]),
                    "source": " | ".join(str(value) for value in candidate.get("sources", [])),
                }
            )
            all_exclusions.append(exclusion)

    all_runs.sort(key=lambda item: (item["study_id"], item["run_accession"]))
    all_samples.sort(key=lambda item: (item["study_id"], item["run_accession"]))
    all_exclusions.sort(key=lambda item: (item["entity_type"], item["study_id"], item["entity_id"]))
    return {
        "studies": studies,
        "runs": all_runs,
        "samples": all_samples,
        "exclusion_log": all_exclusions,
    }
