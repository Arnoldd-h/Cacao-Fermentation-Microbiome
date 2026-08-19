"""Systematic, reproducible discovery of public cacao fermentation BioProjects."""

from __future__ import annotations

import re
from typing import Any

from .schema import DISCOVERY_COLUMNS


def _unique(rows: list[dict[str, str]], field: str) -> str:
    return "; ".join(sorted({row.get(field, "") for row in rows if row.get(field)}))


def _sum_bytes(rows: list[dict[str, str]]) -> int:
    total = 0
    for row in rows:
        for value in row.get("fastq_bytes", "").split(";"):
            if value.strip().isdigit():
                total += int(value.strip())
    return total


def _signals(summary: dict[str, Any], runs: list[dict[str, str]]) -> tuple[bool, bool]:
    text = " ".join(
        [str(summary.get("project_title", "")), str(summary.get("project_description", ""))]
        + [
            " ".join(
                row.get(field, "")
                for field in (
                    "sample_alias",
                    "sample_title",
                    "sample_description",
                    "experiment_title",
                    "library_strategy",
                )
            )
            for row in runs
        ]
    ).lower()
    bacterial_16s = bool(re.search(r"\b16s\b|bacter(?:ia|ial|ium)", text))
    temporal = bool(
        re.search(r"\b\d+(?:\.\d+)?\s*(?:h|hr|hrs|hour|hours|day|days)\b", text)
        or re.search(r"time[- ]?series|longitudinal|succession|during fermentation", text)
    )
    return bacterial_16s, temporal


def preliminary_screen(
    summary: dict[str, Any], runs: list[dict[str, str]], configured: bool
) -> tuple[str, str, bool, bool]:
    """Triage discovery results without promoting them to scientific inclusion."""

    bacterial_16s, temporal = _signals(summary, runs)
    strategies = {row.get("library_strategy", "").upper() for row in runs if row.get("library_strategy")}
    title = str(summary.get("project_title", "")).lower()
    if configured:
        return "configured", "Full screening decision is stored in config/datasets.yaml", bacterial_16s, temporal
    if not runs:
        return "exclude_no_public_runs", "No public ENA read runs were returned", bacterial_16s, temporal
    if strategies and strategies <= {"WGS"}:
        return "exclude_phase_i_wgs", "All discovered runs are WGS, outside Phase I bacterial 16S scope", bacterial_16s, temporal
    if "transcriptome" in title or strategies & {"RNA-SEQ"}:
        return "exclude_phase_i_other", "Transcriptomic study is outside Phase I bacterial 16S scope", bacterial_16s, temporal
    if "AMPLICON" in strategies and bacterial_16s and temporal:
        return "manual_review_priority", "Amplicon, bacterial/16S and temporal signals require full paper-level screening", bacterial_16s, temporal
    if "AMPLICON" in strategies and bacterial_16s:
        return "manual_review_16s", "Amplicon and bacterial/16S signals found; temporal design requires verification", bacterial_16s, temporal
    if "AMPLICON" in strategies:
        return "manual_review_marker", "Amplicon runs found but bacterial 16S marker was not verified automatically", bacterial_16s, temporal
    return "manual_review_technology", "Run technology or marker requires manual verification", bacterial_16s, temporal


def build_discovery_rows(
    search_id: str,
    query: str,
    summaries: dict[str, dict[str, Any]],
    runs_by_project: dict[str, list[dict[str, str]]],
    configured_accessions: set[str],
    retrieved_at_utc: str,
) -> list[dict[str, str]]:
    """Build deterministic discovery rows sorted by accession."""

    output: list[dict[str, str]] = []
    for accession, summary in summaries.items():
        runs = runs_by_project.get(accession, [])
        configured = accession in configured_accessions
        screen, reason, bacterial_16s, temporal = preliminary_screen(summary, runs, configured)
        row = {column: "" for column in DISCOVERY_COLUMNS}
        row.update(
            {
                "search_id": search_id,
                "query": query,
                "bioproject": accession,
                "project_title": str(summary.get("project_title", "")),
                "project_description": str(summary.get("project_description", "")),
                "registration_date": str(summary.get("registration_date", "")),
                "project_data_type": str(summary.get("project_data_type", "")),
                "total_runs": str(len(runs)),
                "library_strategies": _unique(runs, "library_strategy"),
                "library_selections": _unique(runs, "library_selection"),
                "library_layouts": _unique(runs, "library_layout"),
                "sequencing_instruments": _unique(runs, "instrument_model"),
                "countries": _unique(runs, "country"),
                "estimated_fastq_bytes": str(_sum_bytes(runs)),
                "bacterial_16s_signal": str(bacterial_16s).lower(),
                "temporal_signal": str(temporal).lower(),
                "preliminary_screen": screen,
                "screening_reason": reason,
                "already_configured": str(configured).lower(),
                "source": f"https://www.ncbi.nlm.nih.gov/bioproject/{accession} | https://www.ebi.ac.uk/ena/browser/view/{accession}",
                "retrieved_at_utc": retrieved_at_utc,
            }
        )
        output.append(row)
    output.sort(key=lambda item: item["bioproject"])
    return output
