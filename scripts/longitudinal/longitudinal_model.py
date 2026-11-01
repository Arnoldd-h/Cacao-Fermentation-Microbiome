"""Deterministic CLR trajectories and independent-batch slope summaries.

All analytical settings come from the registered longitudinal configuration.
No time permutation, ASV selection, interpolation or sample-level inference.
"""

import itertools
import math
from collections import defaultdict

from cacao_inventory.analysis_artifacts import unique_rows

IDENTITY = ["study_id", "sample_id", "run_accession", "fermentation_batch", "fermentation_hours",
            "relative_time", "fermentation_duration_hours", "sampling_stratum"]
OBSERVATION_FIELDS = ["scenario_id", "study_id", "fermentation_batch", "fermentation_hours", "relative_time",
                      "time_value", "time_unit", "pseudocount", "n_strata", "source_sample_ids", "summed_source_reads", "batch_time_index"]
MEMBERSHIP_FIELDS = ["scenario_id", *IDENTITY, "included", "reason"]
MODEL_FIELDS = ["scenario_id", "study_id", "fermentation_batch", "n_timepoints", "n_samples", "time_unit",
                "time_min", "time_max", "slope_norm", "residual_sum_squares", "total_sum_squares", "r_squared", "response_status"]
GLOBAL_FIELDS = ["scenario_id", "study_id", "time_unit", "n_batches", "n_observations", "n_samples",
                 "sign_patterns", "minimum_attainable_p", "norm_mean_slope", "mean_batch_slope_norm",
                 "between_batch_rms", "exact_p_value", "fdr_by_q_value", "reject_fdr"]
DIRECTION_FIELDS = ["scenario_id", "batch_1", "batch_2", "cosine", "status"]
LOO_FIELDS = ["scenario_id", "omitted_batch", "remaining_batches", "norm_mean_slope", "cosine_with_all_batches", "status"]
SIGN_FIELDS = ["scenario_id", "batch_signs", "statistic", "at_least_observed"]
PLOT_FIELDS = ["study_id", "fermentation_batch", "fermentation_hours", "relative_time", "batch_time_index", "n_strata", "distance_from_baseline"]


def vector_norm(vector):
    return math.sqrt(math.fsum(value * value for value in vector))


def mean_vector(vectors):
    if not vectors or not vectors[0] or any(len(row) != len(vectors[0]) for row in vectors):
        raise ValueError("Cannot average empty or incompatible vectors")
    return [math.fsum(values) / len(vectors) for values in zip(*vectors)]


def cosine(left, right):
    denominator = vector_norm(left) * vector_norm(right)
    if denominator == 0:
        return "", "zero_slope"
    return max(-1.0, min(1.0, math.fsum(a * b for a, b in zip(left, right)) / denominator)), "defined"


def clr_vector(counts, pseudocount):
    if not counts or any(type(value) is not int or value < 0 for value in counts) or sum(counts) <= 0:
        raise ValueError("CLR requires nonnegative integer counts and a nonempty library")
    if type(pseudocount) not in (int, float) or not math.isfinite(pseudocount) or pseudocount <= 0:
        raise ValueError("Invalid pseudocount")
    logs = [math.log(value + pseudocount) for value in counts]
    center = math.fsum(logs) / len(logs)
    return [value - center for value in logs]


def validate_cohort(counts, features, metadata, config):
    cohort = config["cohort"]
    samples = unique_rows(counts, ("study_id", "sample_id"))
    meta = unique_rows(metadata, ("study_id", "sample_id"))
    if samples.keys() != meta.keys() or len(samples) != cohort["expected_samples"] or len(features) < 2:
        raise ValueError("Longitudinal sample/feature scope mismatch")
    if {row["study_id"] for row in counts} != {cohort["study_id"]}:
        raise ValueError("Mixed or unexpected study; ASVs cannot be pooled")
    if len({row["run_accession"] for row in metadata}) != len(metadata):
        raise ValueError("Repeated run accession")
    groups = defaultdict(list)
    identities = set()
    tolerance = config["method"]["validation_tolerance"]
    for key, row in meta.items():
        if any(not row[field] for field in IDENTITY):
            raise ValueError("Missing longitudinal identity or time metadata")
        hours, relative, duration = (float(row[field]) for field in ("fermentation_hours", "relative_time", "fermentation_duration_hours"))
        if not all(math.isfinite(value) for value in (hours, relative, duration)) or hours < 0 or duration <= 0 or not 0 <= relative <= 1:
            raise ValueError("Invalid time or duration")
        if not math.isclose(relative, hours / duration, rel_tol=tolerance, abs_tol=tolerance):
            raise ValueError("Relative time disagrees with source hours/duration")
        identity = (row["fermentation_batch"], hours, row["sampling_stratum"])
        if identity in identities or row["sampling_stratum"] not in cohort["expected_strata"]:
            raise ValueError("Duplicate batch/time/stratum or unexpected stratum")
        identities.add(identity)
        groups[identity[:2]].append(row)
        if sum(int(samples[key][feature]) for feature in features) <= 0:
            raise ValueError("Empty longitudinal library")
    batches = {batch for batch, _ in groups}
    if len(batches) != cohort["expected_batches"] or len(groups) != cohort["expected_batch_times"]:
        raise ValueError("Batch/time scope mismatch")
    for batch in batches:
        rows = [row for (current, _), values in groups.items() if current == batch for row in values]
        durations = {float(row["fermentation_duration_hours"]) for row in rows}
        if len(durations) != 1 or not math.isclose(next(iter(durations)), max(float(row["fermentation_hours"]) for row in rows), rel_tol=tolerance, abs_tol=tolerance):
            raise ValueError("Duration is not the registered last observed time")
    for values in groups.values():
        if len({float(row["relative_time"]) for row in values}) != 1:
            raise ValueError("Subsamples disagree on relative time")
    return samples, groups


def aggregate_observations(counts, features, metadata, config, scenario):
    samples, groups = validate_cohort(counts, features, metadata, config)
    full_strata = set(config["cohort"]["expected_strata"])
    indexes = {}
    for batch in sorted({key[0] for key in groups}):
        indexes.update({(batch, hour): index + 1 for index, hour in enumerate(sorted(hour for current, hour in groups if current == batch))})
    observations, membership, vectors = [], [], {}
    for (batch, hours), values in sorted(groups.items()):
        values = sorted(values, key=lambda row: (row["sampling_stratum"], row["sample_id"]))
        complete = {row["sampling_stratum"] for row in values} == full_strata
        selected = []
        for row in values:
            include, reason = True, "included"
            if scenario["strata_policy"] == "complete" and not complete:
                include, reason = False, "incomplete_strata_observation"
            elif scenario["strata_policy"] == "single" and row["sampling_stratum"] != scenario["stratum"]:
                include, reason = False, "other_stratum_sensitivity"
            membership.append({"scenario_id": scenario["id"], **{field: row[field] for field in IDENTITY},
                               "included": "TRUE" if include else "FALSE", "reason": reason})
            if include:
                selected.append(row)
        if not selected:
            continue
        sample_counts = [[int(samples[row["study_id"], row["sample_id"]][feature]) for feature in features] for row in selected]
        vector = mean_vector([clr_vector(row, scenario["pseudocount"]) for row in sample_counts])
        vectors[batch, hours] = vector
        relative = float(selected[0]["relative_time"])
        time = relative if scenario["time_field"] == "relative_time" else hours
        observations.append({"scenario_id": scenario["id"], "study_id": config["cohort"]["study_id"],
            "fermentation_batch": batch, "fermentation_hours": hours, "relative_time": relative,
            "time_value": time / scenario["time_scale"], "time_unit": scenario["time_unit"], "pseudocount": scenario["pseudocount"],
            "n_strata": len(selected), "source_sample_ids": ";".join(row["sample_id"] for row in selected),
            "summed_source_reads": sum(sum(row) for row in sample_counts), "batch_time_index": indexes[batch, hours]})
    for batch in sorted({key[0] for key in groups}):
        series = [row for row in observations if row["fermentation_batch"] == batch]
        if len(series) < config["cohort"]["minimum_timepoints_per_batch"] or not any(row["fermentation_hours"] == 0 for row in series):
            raise ValueError("Every scenario requires all batches, baseline and enough timepoints")
    return observations, vectors, membership


def fit_multivariate_ols(times, vectors):
    if len(times) != len(vectors) or len(times) < 3 or len(set(times)) != len(times):
        raise ValueError("OLS requires at least three distinct, matched times")
    average_x = math.fsum(times) / len(times)
    centered_x = [value - average_x for value in times]
    denominator = math.fsum(value * value for value in centered_x)
    if denominator <= 0:
        raise ValueError("Temporal model is rank deficient")
    means = mean_vector(vectors)
    slopes = [0.0 if len(set(values)) == 1 else math.fsum(dx * (y - mean) for dx, y in zip(centered_x, values)) / denominator
              for values, mean in zip(zip(*vectors), means)]
    intercepts = [mean - average_x * slope for mean, slope in zip(means, slopes)]
    residual = math.fsum((value - (intercept + time * slope))**2 for time, row in zip(times, vectors)
                         for value, intercept, slope in zip(row, intercepts, slopes))
    total = math.fsum((value - mean)**2 for row in vectors for value, mean in zip(row, means))
    constant = all(row == vectors[0] for row in vectors)
    return slopes, intercepts, residual, total, "" if constant else 1 - residual / total


def exact_sign_test(slopes, tolerance, maximum_patterns):
    size = 2 ** len(slopes)
    if size > maximum_patterns:
        raise ValueError("Sign enumeration exceeds the configured limit")
    observed = vector_norm(mean_vector(slopes))
    threshold = observed - tolerance * max(1, observed)
    patterns = []
    for signs in itertools.product((-1, 1), repeat=len(slopes)):
        statistic = vector_norm([math.fsum(sign * row[j] for sign, row in zip(signs, slopes)) / len(slopes) for j in range(len(slopes[0]))])
        patterns.append({"batch_signs": "".join("+" if sign == 1 else "-" for sign in signs),
                         "statistic": statistic, "at_least_observed": "TRUE" if statistic >= threshold else "FALSE"})
    p_value = sum(row["at_least_observed"] == "TRUE" for row in patterns) / size
    return observed, p_value, patterns


def adjust_by(p_values):
    if not p_values or any(not math.isfinite(value) or not 0 <= value <= 1 for value in p_values):
        raise ValueError("Invalid p-values")
    count = len(p_values)
    harmonic = math.fsum(1 / rank for rank in range(1, count + 1))
    order = sorted(range(count), key=lambda index: p_values[index])
    result, previous = [0.0] * count, 1.0
    for position in range(count - 1, -1, -1):
        index = order[position]
        previous = min(previous, count * harmonic * p_values[index] / (position + 1))
        result[index] = previous
    return result


def build_tables(counts, features, metadata, config):
    tables = {name: [] for name in ("observations.tsv", "sample_membership.tsv", "batch_models.tsv", "slopes.tsv", "intercepts.tsv",
                                   "global_tests.tsv", "direction_cosines.tsv", "leave_one_batch_out.tsv", "sign_statistics.tsv", "trajectory_plot_data.tsv")}
    study = config["cohort"]["study_id"]
    for scenario in config["scenarios"]:
        observations, vectors, membership = aggregate_observations(counts, features, metadata, config, scenario)
        tables["observations.tsv"].extend(observations)
        tables["sample_membership.tsv"].extend(membership)
        batches = sorted({row["fermentation_batch"] for row in observations})
        slopes = []
        for batch in batches:
            rows = [row for row in observations if row["fermentation_batch"] == batch]
            response = [vectors[batch, row["fermentation_hours"]] for row in rows]
            slope, intercept, residual, total, r_squared = fit_multivariate_ols([row["time_value"] for row in rows], response)
            slopes.append(slope)
            identity = {"scenario_id": scenario["id"], "study_id": study, "fermentation_batch": batch}
            tables["slopes.tsv"].append({**identity, **dict(zip(features, slope))})
            tables["intercepts.tsv"].append({**identity, **dict(zip(features, intercept))})
            tables["batch_models.tsv"].append({**identity, "n_timepoints": len(rows), "n_samples": sum(row["n_strata"] for row in rows),
                "time_unit": scenario["time_unit"], "time_min": min(row["time_value"] for row in rows), "time_max": max(row["time_value"] for row in rows),
                "slope_norm": vector_norm(slope), "residual_sum_squares": residual, "total_sum_squares": total,
                "r_squared": r_squared, "response_status": "constant" if r_squared == "" else "variable"})
            if scenario["id"] == config["primary_scenario"]:
                baseline = vectors[batch, 0.0]
                tables["trajectory_plot_data.tsv"].extend({**{field: row[field] for field in PLOT_FIELDS if field != "distance_from_baseline"},
                    "distance_from_baseline": vector_norm([a - b for a, b in zip(vectors[batch, row["fermentation_hours"]], baseline)])} for row in rows)
        average = mean_vector(slopes)
        observed, p_value, patterns = exact_sign_test(slopes, config["method"]["comparison_tolerance"], config["method"]["maximum_sign_patterns"])
        tables["sign_statistics.tsv"].extend({"scenario_id": scenario["id"], **row} for row in patterns)
        tables["global_tests.tsv"].append({"scenario_id": scenario["id"], "study_id": study, "time_unit": scenario["time_unit"],
            "n_batches": len(batches), "n_observations": len(observations), "n_samples": sum(row["n_strata"] for row in observations),
            "sign_patterns": len(patterns), "minimum_attainable_p": 2 / len(patterns), "norm_mean_slope": observed,
            "mean_batch_slope_norm": math.fsum(map(vector_norm, slopes)) / len(slopes),
            "between_batch_rms": math.sqrt(math.fsum((value - mean)**2 for row in slopes for value, mean in zip(row, average)) / len(slopes)),
            "exact_p_value": p_value})
        for left, right in itertools.combinations(range(len(batches)), 2):
            value, status = cosine(slopes[left], slopes[right])
            tables["direction_cosines.tsv"].append({"scenario_id": scenario["id"], "batch_1": batches[left], "batch_2": batches[right], "cosine": value, "status": status})
        for index, batch in enumerate(batches):
            remaining = mean_vector([row for position, row in enumerate(slopes) if position != index])
            value, status = cosine(remaining, average)
            tables["leave_one_batch_out.tsv"].append({"scenario_id": scenario["id"], "omitted_batch": batch,
                "remaining_batches": len(batches) - 1, "norm_mean_slope": vector_norm(remaining), "cosine_with_all_batches": value, "status": status})
    for row, q_value in zip(tables["global_tests.tsv"], adjust_by([row["exact_p_value"] for row in tables["global_tests.tsv"]])):
        row.update(fdr_by_q_value=q_value, reject_fdr="TRUE" if q_value <= config["method"]["alpha"] else "FALSE")
    return tables


def table_columns(features):
    coefficients = ["scenario_id", "study_id", "fermentation_batch", *features]
    return dict(zip(("observations.tsv", "sample_membership.tsv", "batch_models.tsv", "slopes.tsv", "intercepts.tsv", "global_tests.tsv",
                     "direction_cosines.tsv", "leave_one_batch_out.tsv", "sign_statistics.tsv", "trajectory_plot_data.tsv"),
                    (OBSERVATION_FIELDS, MEMBERSHIP_FIELDS, MODEL_FIELDS, coefficients, coefficients, GLOBAL_FIELDS,
                     DIRECTION_FIELDS, LOO_FIELDS, SIGN_FIELDS, PLOT_FIELDS)))
