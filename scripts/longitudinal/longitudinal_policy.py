"""Validate explicit within-study temporal methods before reading abundances."""

import math
import re

from cacao_inventory.taxonomy_reference import repository_path


def positive_number(value):
    return type(value) in (int, float) and math.isfinite(value) and value > 0


def validate_longitudinal_policy(config, root):
    source = repository_path(root, config["source_diversity_config"])
    output = repository_path(root, config["output_dir"])
    if not source.is_relative_to(root.resolve() / "config") or not output.is_relative_to(root.resolve() / "results/longitudinal"):
        raise ValueError("Longitudinal config/output paths outside their registered roots")
    cohort, method = config["cohort"], config["method"]
    for field in ("expected_samples", "expected_batches", "expected_batch_times", "minimum_timepoints_per_batch"):
        if type(cohort[field]) is not int or cohort[field] < 1:
            raise ValueError(f"Invalid cohort {field}")
    if cohort["expected_batches"] < 2 or cohort["minimum_timepoints_per_batch"] < 3 or not cohort["study_id"]:
        raise ValueError("Independent batches and at least three timepoints are required")
    strata = cohort["expected_strata"]
    if len(strata) < 2 or len(set(strata)) != len(strata) or any(not isinstance(s, str) or not s for s in strata):
        raise ValueError("Invalid expected sampling strata")
    if cohort["duration_basis"] != "last_observed_time_in_registered_metadata":
        raise ValueError("Unsupported duration policy")
    expected = {
        "scope": "exploratory_within_study", "independent_unit": "fermentation_batch",
        "feature_level": "within_study_asv", "zero_policy": "add_pseudocount_to_all_counts",
        "aggregation": "equal_weight_mean_clr_within_batch_time",
        "model": "multivariate_ordinary_least_squares_with_intercept_per_batch",
        "combination": "equal_weight_mean_batch_slope",
        "global_test": "exhaustive_batch_sign_flip_norm_mean_slope",
        "null_assumptions": ["independent_batches", "centrally_symmetric_batch_slope_vectors_about_zero"],
        "fdr_method": "BY", "fdr_family": "all_configured_global_scenarios",
        "confidence_intervals": False, "uncertainty": "batch_effects_and_leave_one_batch_out",
        "pilot_inference": False, "random_sampling": False}
    if any(method.get(key) != value or (isinstance(value, bool) and type(method.get(key)) is not bool) for key, value in expected.items()):
        raise ValueError("Unsupported longitudinal method or independence policy")
    if type(method["random_seed"]) is not int or method["random_seed"] < 0 or type(method["maximum_sign_patterns"]) is not int:
        raise ValueError("Invalid seed or enumeration limit")
    if 2 ** cohort["expected_batches"] > method["maximum_sign_patterns"]:
        raise ValueError("Exact enumeration exceeds registered limit; no random approximation")
    if not positive_number(method["alpha"]) or method["alpha"] >= 1:
        raise ValueError("Invalid alpha")
    for field in ("comparison_tolerance", "validation_tolerance"):
        if not positive_number(method[field]) or method[field] > 1e-6:
            raise ValueError(f"Invalid {field}")
    scenarios = config["scenarios"]
    ids = [row["id"] for row in scenarios]
    if not scenarios or len(set(ids)) != len(ids) or config["primary_scenario"] not in ids:
        raise ValueError("Missing/duplicate scenarios or primary scenario")
    for row in scenarios:
        if not re.fullmatch(r"[a-z][a-z0-9_]*", row["id"]):
            raise ValueError("Invalid scenario identifier")
        if row["time_field"] not in ("relative_time", "fermentation_hours") or not row["time_unit"]:
            raise ValueError("Unknown time field or missing time unit")
        if not positive_number(row["time_scale"]) or not positive_number(row["pseudocount"]):
            raise ValueError("Time scale and pseudocount must be finite and positive")
        if row["strata_policy"] not in ("complete", "available", "single"):
            raise ValueError("Unknown strata aggregation policy")
        if ((row["strata_policy"] == "single" and row.get("stratum") not in strata) or
                (row["strata_policy"] != "single" and "stratum" in row)):
            raise ValueError("Single-stratum policy requires a registered stratum")
    figures = config["figures"]
    if any(not positive_number(figures[key]) for key in ("width_inches", "height_inches")) or type(figures["png_dpi"]) is not int or figures["png_dpi"] < 300:
        raise ValueError("Invalid publication figure dimensions/DPI")
    if not re.fullmatch(r"#[0-9a-fA-F]{6}", figures["color"]):
        raise ValueError("Invalid figure color")
    return source, output
