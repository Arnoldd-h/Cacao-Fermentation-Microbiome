"""Validate study membership and promotion criteria after paired DADA2.

This module reads small tables only. It does not replace FASTQ/ASV validation,
provenance checks or a scientific interpretation of the resulting sequences.
"""

import csv
from decimal import Decimal, InvalidOperation
from pathlib import Path
import re


NUMERIC_METADATA = frozenset({
    "pilot_order", "fermentation_hours", "relative_time",
    "fermentation_duration_hours", "estimated_bytes_total",
})
IDENTITY_FIELDS = ("study_id", "sample_id", "run_accession")
TRACKING_STAGES = ("input", "filtered", "denoised_f", "denoised_r", "merged", "nonchim")


def _require(row: dict, fields, context: str):
    missing = set(fields) - set(row)
    if missing:
        raise ValueError(f"Missing fields in {context}: {sorted(missing)}")


def _integer(value, context: str) -> int:
    if isinstance(value, bool) or not re.fullmatch(r"[0-9]+", str(value)):
        raise ValueError(f"Expected a nonnegative integer for {context}: {value!r}")
    return int(value)


def _number(value, context: str) -> Decimal:
    try:
        number = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError(f"Invalid numeric metadata for {context}: {value!r}") from exc
    if not number.is_finite():
        raise ValueError(f"Nonfinite numeric metadata for {context}")
    return number


def _same_metadata_value(field: str, observed, expected) -> bool:
    if field in NUMERIC_METADATA and observed != "" and expected != "":
        return _number(observed, field) == _number(expected, field)
    return observed == expected


def _unique(rows: list[dict], field: str, context: str) -> dict:
    result = {}
    if not rows:
        raise ValueError(f"Empty {context}")
    for row in rows:
        _require(row, [field], context)
        key = row[field]
        if not isinstance(key, str) or not key.strip() or key in result:
            raise ValueError(f"Missing or duplicate {field} in {context}: {key!r}")
        result[key] = row
    return result


def _manifest_index(manifest: list[dict]) -> dict:
    samples = _unique(manifest, "sample_id", "manifest")
    _unique(manifest, "run_accession", "manifest")
    for row in manifest:
        _require(row, IDENTITY_FIELDS, "manifest")
        if not row["study_id"]:
            raise ValueError("Empty study_id in manifest")
    if len({row["study_id"] for row in manifest}) != 1:
        raise ValueError("Manifest mixes studies")
    return samples


def _read_table(directory: Path, name: str) -> list[dict]:
    with (directory / name).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames
        if not columns or len(columns) != len(set(columns)) or any(not field for field in columns):
            raise ValueError(f"Missing or duplicate table columns: {name}")
        rows = list(reader)
    if not rows or any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Empty or malformed table: {name}")
    return rows


def quality_rows(manifest: list[dict], terminal_rows: list[dict]) -> list[dict]:
    """Return the R runner's paired QC-count rows in manifest, then R1/R2 order.

    Terminal rows use ``direction`` and ``pairs``. Membership and identities
    must match the manifest exactly, with one row per sample and mate.
    """
    samples = _manifest_index(manifest)
    observed = {}
    for row in terminal_rows:
        _require(row, (*IDENTITY_FIELDS, "direction", "pairs"), "terminal QC")
        sample = row["sample_id"]
        key = (sample, row["direction"])
        if sample not in samples or row["direction"] not in ("R1", "R2") or key in observed:
            raise ValueError("Terminal QC has duplicate, extra or invalid sample/mate membership")
        if any(row[field] != samples[sample][field] for field in IDENTITY_FIELDS):
            raise ValueError(f"Terminal QC identity disagrees with manifest: {sample}")
        observed[key] = _integer(row["pairs"], f"terminal QC pairs for {sample}")
    expected = {(sample, direction) for sample in samples for direction in ("R1", "R2")}
    if set(observed) != expected:
        raise ValueError("Terminal QC sample/mate membership is incomplete")
    result = []
    for row in manifest:
        sample = row["sample_id"]
        if observed[sample, "R1"] != observed[sample, "R2"]:
            raise ValueError(f"Terminal QC paired counts disagree: {sample}")
        for direction in ("R1", "R2"):
            result.append({
                **{field: row[field] for field in IDENTITY_FIELDS},
                "read_direction": direction, "total_sequences": observed[sample, direction],
            })
    return result


def _validate_metadata(directory: Path, manifest: list[dict], samples: dict):
    metadata = _unique(_read_table(directory, "sample_metadata.tsv"), "sample_id", "sample metadata")
    if set(metadata) != set(samples):
        raise ValueError("Sample metadata membership differs from manifest")
    columns = set(manifest[0])
    for sample, expected in samples.items():
        observed = metadata[sample]
        if set(expected) != columns or set(observed) != columns:
            raise ValueError("Sample metadata must preserve every manifest field exactly")
        for field in columns:
            if not _same_metadata_value(field, observed[field], expected[field]):
                raise ValueError(f"Sample metadata altered {field}: {sample}")


def _cohort_design(samples: dict, policy: dict) -> tuple[int, int, int]:
    batches, batch_times = set(), set()
    unknown_strata = 0
    for sample, row in samples.items():
        _require(row, ("fermentation_batch", "fermentation_hours", "sampling_stratum"), "manifest")
        if not row["fermentation_batch"]:
            raise ValueError(f"Missing fermentation batch: {sample}")
        hours = _number(row["fermentation_hours"], "fermentation_hours")
        if hours < 0:
            raise ValueError(f"Negative fermentation time: {sample}")
        batches.add(row["fermentation_batch"])
        batch_times.add((row["fermentation_batch"], hours))
        unknown_strata += row["sampling_stratum"] == ""
    for field, observed in (("expected_batches", len(batches)), ("expected_batch_times", len(batch_times))):
        _require(policy, [field], "cohort policy")
        if _integer(policy[field], field) != observed:
            raise ValueError(f"Cohort design disagrees with {field}")
    if "expected_runs" in policy and _integer(policy["expected_runs"], "expected_runs") != len(samples):
        raise ValueError("Cohort design disagrees with expected_runs")
    return len(batches), len(batch_times), unknown_strata


def _tracking_counts(directory: Path, samples: dict, inputs: dict) -> tuple[dict, dict]:
    tracking = _unique(_read_table(directory, "read_tracking.tsv"), "sample_id", "read tracking")
    if set(tracking) != set(samples):
        raise ValueError("Read tracking membership differs from manifest")
    counts, totals = {}, {stage: 0 for stage in TRACKING_STAGES}
    for sample, row in tracking.items():
        _require(row, (*IDENTITY_FIELDS, *TRACKING_STAGES), "read tracking")
        for field in set(row) & set(samples[sample]):
            if not _same_metadata_value(field, row[field], samples[sample][field]):
                raise ValueError(f"Tracking metadata altered {field}: {sample}")
        count = {stage: _integer(row[stage], f"{sample} {stage}") for stage in TRACKING_STAGES}
        if count["input"] != inputs[sample]:
            raise ValueError(f"DADA2 input differs from terminal QC: {sample}")
        if "raw_input" in row and _integer(row["raw_input"], "raw_input") != inputs[sample]:
            raise ValueError(f"DADA2 raw_input differs from terminal QC: {sample}")
        if not (0 < count["nonchim"] <= count["merged"] <= min(count["denoised_f"], count["denoised_r"])
                <= max(count["denoised_f"], count["denoised_r"]) <= count["filtered"] <= count["input"]):
            raise ValueError(f"Read conservation or positive nonchim requirement failed: {sample}")
        counts[sample] = count
        for stage in TRACKING_STAGES:
            totals[stage] += count[stage]
    return counts, totals


def _validate_merges(directory: Path, samples: dict, counts: dict):
    merges = _unique(_read_table(directory, "merge_diagnostics.tsv"), "sample_id", "merge diagnostics")
    if set(merges) != set(samples):
        raise ValueError("Merge diagnostic membership differs from manifest")
    for sample, row in merges.items():
        _require(row, ("study_id", "accepted_pairs", "rejected_pairs"), "merge diagnostics")
        if row["study_id"] != samples[sample]["study_id"]:
            raise ValueError(f"Merge diagnostic study differs from manifest: {sample}")
        if "run_accession" in row and row["run_accession"] != samples[sample]["run_accession"]:
            raise ValueError(f"Merge diagnostic run differs from manifest: {sample}")
        accepted = _integer(row["accepted_pairs"], "accepted_pairs")
        _integer(row["rejected_pairs"], "rejected_pairs")
        if accepted != counts[sample]["merged"]:
            raise ValueError(f"Accepted merge pairs differ from tracking: {sample}")


def _boolean(value: str, context: str) -> bool:
    if value not in ("TRUE", "FALSE"):
        raise ValueError(f"Expected TRUE/FALSE in {context}")
    return value == "TRUE"


def _error_diagnostics(directory: Path) -> dict:
    errors = _unique(_read_table(directory, "error_learning.tsv"), "read_direction", "error learning")
    if set(errors) != {"F", "R"}:
        raise ValueError("Error learning requires exactly F and R diagnostics")
    result = {}
    for direction, row in errors.items():
        _require(row, ("self_consistency_rounds", "maximum_rounds", "repeated_error_matrix",
                       "reached_iteration_limit", "converged_before_limit"), "error learning")
        rounds = _integer(row["self_consistency_rounds"], "self_consistency_rounds")
        maximum = _integer(row["maximum_rounds"], "maximum_rounds")
        repeated = _boolean(row["repeated_error_matrix"], "repeated_error_matrix")
        reached = _boolean(row["reached_iteration_limit"], "reached_iteration_limit")
        converged = _boolean(row["converged_before_limit"], "converged_before_limit")
        if not 1 <= rounds <= maximum or reached != (rounds == maximum) or converged != (repeated and rounds < maximum):
            raise ValueError(f"Inconsistent error-learning rounds or convergence flags: {direction}")
        if "transition_base_observations" in row:
            _integer(row["transition_base_observations"], "transition_base_observations")
        result[direction] = {
            "self_consistency_rounds": rounds, "maximum_rounds": maximum,
            "repeated_error_matrix": repeated, "reached_iteration_limit": reached,
            "converged_before_limit": converged,
        }
    return result


def validate_cohort(directory: Path, manifest: list[dict], terminal_rows: list[dict],
                    terminal_summary: dict, policy: dict) -> dict:
    """Check complete cohort preservation and return a taxonomy promotion gate.

    Inputs are the DADA2 output directory, the full source manifest, terminal QC
    run-summary rows/summary and a policy with expected_batches/batch_times.
    Structural/count/metadata inconsistencies raise ValueError. Valid outputs
    whose error model did not converge before its limit return needs_review.
    """
    directory = Path(directory)
    samples = _manifest_index(manifest)
    rows = quality_rows(manifest, terminal_rows)
    inputs = {row["sample_id"]: row["total_sequences"] for row in rows if row["read_direction"] == "R1"}
    study = manifest[0]["study_id"]
    _require(terminal_summary, ("study_id", "runs", "pairs", "pairs_below_min_len"), "terminal summary")
    if terminal_summary["study_id"] != study or _integer(terminal_summary["runs"], "terminal runs") != len(samples):
        raise ValueError("Terminal summary cohort differs from manifest")
    total_pairs = _integer(terminal_summary["pairs"], "terminal pairs")
    below = _integer(terminal_summary["pairs_below_min_len"], "terminal pairs_below_min_len")
    if total_pairs != sum(inputs.values()) or below > total_pairs:
        raise ValueError("Terminal summary pair counts disagree with cohort")
    _validate_metadata(directory, manifest, samples)
    batches, batch_times, unknown_strata = _cohort_design(samples, policy)
    counts, totals = _tracking_counts(directory, samples, inputs)
    eligible = total_pairs - below
    if totals["filtered"] > eligible:
        raise ValueError("DADA2 filtered pairs exceed terminal minimum-length eligibility")
    _validate_merges(directory, samples, counts)
    errors = _error_diagnostics(directory)
    review_directions = [direction for direction in ("F", "R") if not errors[direction]["converged_before_limit"]]
    return {
        "status": "valid", "study_id": study, "samples": len(samples),
        "batches": batches, "batch_times": batch_times,
        "unknown_sampling_stratum_samples": unknown_strata,
        "input_pairs": totals["input"], "filtered_pairs": totals["filtered"],
        "merged_pairs": totals["merged"], "nonchim_pairs": totals["nonchim"],
        "pairs_at_least_min_len": eligible,
        "minimum_nonchim_pairs_per_sample": min(row["nonchim"] for row in counts.values()),
        "retention": {
            "filtered_fraction_of_input": totals["filtered"] / total_pairs,
            "merged_fraction_of_filtered": totals["merged"] / totals["filtered"],
            "nonchim_fraction_of_merged": totals["nonchim"] / totals["merged"],
            "nonchim_fraction_of_input": totals["nonchim"] / total_pairs,
        },
        "error_learning": errors,
        "promotion_status": "needs_review" if review_directions else "ready_for_taxonomy",
        "promotion_reason": ("Error models did not converge before the iteration limit: " + ", ".join(review_directions)
                             if review_directions else "Complete cohort and both error models converged before their limits"),
    }
