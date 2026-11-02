#!/usr/bin/env python3
"""Execute or validate independent DADA2 after validated terminal QC.

Inputs: execution policy, original processing profile/manifest, validated QC
and paired terminal FASTQ. Outputs: DADA2 tables/figures, separate SHA-256
audit and ignored filtered FASTQ/RDS/log. Fails on stale input hashes, changed
cohort/counts/settings/software or unavailable inputs. Nonconvergence preserves
technically valid results with needs_review, without automatic parameter tuning.
Existing successful executions are validated, never silently overwritten.
"""
import argparse
import json
import logging
import math
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/dada2"))
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, read_table, validate_artifacts, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.pilot_workflow import load_pilot_manifest, processing_configuration
from cacao_inventory.provenance import collect_versions, file_records
from cacao_inventory.taxonomy_reference import repository_path
from validate_outputs import validate_run
from cohort import quality_rows, validate_cohort

HELPERS = ("analysis_artifacts", "config", "download", "io", "pilot_workflow", "provenance", "schema", "taxonomy_reference")
DADA_PRODUCTS = ["SUCCESS", "asv_counts.tsv", "asv_sequences.tsv", "asv_sequences.fasta",
    "read_tracking.tsv", "summary.tsv", "sample_metadata.tsv", "sequence_length_distribution.tsv",
    "error_learning.tsv", "merge_diagnostics.tsv", "input_checksums.tsv", "config_snapshot.yaml",
    "software_versions.tsv", "output_checksums.tsv", "run_provenance.yaml", "warnings.txt",
    "session_info.txt", "validation.json",
    *[f"error_model_{mate}.{ext}" for mate in ("F", "R") for ext in ("pdf", "svg", "png")]]
AUDIT_PRODUCTS = ["error_diagnostics.tsv", "cohort_validation.json", "execution_log.txt"]
TERMINAL_PRODUCTS = ["run_summary.tsv", "read_length_counts.tsv", "quality_coverage.tsv", "summary.json",
                     "config_snapshot.yaml", "input_checksums.json"]


def load_policy(config_path, root=ROOT):
    policy = load_json_yaml(config_path)
    required = {"schema_version": "1.0.0", "method": "independent_study_registered_dada2",
        "processing_ceiling": "dada2", "require_convergence_before_limit": True,
        "require_last_error_matrix_match": True, "sample_exclusions": False,
        "parameter_overrides": False,
        "residual_primer_policy": "diagnose_separately_without_posthoc_exclusions"}
    if any(type(policy.get(key)) is not type(value) or policy[key] != value for key, value in required.items()):
        raise ValueError("Unsupported independent DADA2 policy")
    terminal = load_json_yaml(repository_path(root, policy["source_terminal_policy"]))
    scope = load_json_yaml(repository_path(root, terminal["source_scope"]))
    project = load_json_yaml(repository_path(root, scope["project_config"]))
    manifest = load_pilot_manifest(repository_path(root, scope["manifest"]))
    processing = processing_configuration(project, manifest[0])
    if (policy["expected_runs"] != scope["expected_runs"] or len(manifest) != policy["expected_runs"]
            or policy["expected_batches"] != scope["expected_batches"]
            or any(row["study_id"] != scope["study_id"] or row["bioproject"] != scope["bioproject"] for row in manifest)):
        raise ValueError("Independent DADA2 requires the complete configured cohort")
    for key in ("expected_runs", "expected_batches", "expected_batch_times", "threads", "memory_mb"):
        if type(policy[key]) is not int or policy[key] < 1:
            raise ValueError(f"Invalid positive integer policy field: {key}")
    if (policy["threads"] > processing["dada2"]["threads"]
            or policy["memory_mb"] > scope["resources"]["memory_budget_mb"]
            or not isinstance(policy["minimum_free_space_factor"], (int, float))
            or not math.isfinite(policy["minimum_free_space_factor"])
            or policy["minimum_free_space_factor"] < 3):
        raise ValueError("Invalid resource limits")
    paths = []
    for key, parent in (("output_dir", "results/dada2"), ("audit_dir", "results/dada2"),
                        ("intermediate_dir", "results/intermediate/dada2")):
        path = repository_path(root, policy[key])
        base = (root / parent / scope["study_id"]).resolve()
        if not path.is_relative_to(base) or path == base:
            raise ValueError("Invalid DADA2 output namespace")
        paths.append(path)
    protected = [repository_path(root, scope[key]) for key in ("dada2_dir", "dada2_intermediate_dir")]
    for index, path in enumerate(paths):
        if any(path == other or path in other.parents or other in path.parents
               for other in [*paths[index + 1:], *protected]):
            raise ValueError("DADA2 outputs must not overlap each other or baseline outputs")
    return policy, terminal, scope, manifest, processing


def inputs(config_path, policy, terminal, scope, manifest):
    qc = ROOT / terminal["output_dir"]
    return [config_path, ROOT / policy["source_terminal_policy"], ROOT / terminal["source_scope"],
        ROOT / scope["project_config"], ROOT / scope["manifest"],
        ROOT / "environment/conda-linux-64.lock", ROOT / "workflow/study_dada2.Snakefile",
        *list((ROOT / "scripts/study_dada2").glob("*.py")), ROOT / "scripts/study_dada2/error_diagnostics.R",
        *[ROOT / "scripts/dada2" / name for name in ("run_pilot.R", "helpers.R", "validate_outputs.py")],
        *[ROOT / "python/cacao_inventory" / f"{name}.py" for name in HELPERS],
        *[qc / name for name in ("SUCCESS", "validation.json", "provenance.json", *TERMINAL_PRODUCTS)],
        ROOT / scope["qc_dir"] / "raw_read_quality.tsv",
        ROOT / scope["qc_dir"] / "provenance.json",
        *[ROOT / terminal["interim_dir"] / row["run_accession"] / f'{row["run_accession"]}_{mate}.fastq.gz'
          for row in manifest for mate in (1, 2)]]


def dada_command(policy, terminal, scope):
    return ["Rscript", "scripts/dada2/run_pilot.R", "--config", scope["project_config"],
        "--manifest", scope["manifest"], "--input-dir", terminal["interim_dir"],
        "--output-dir", policy["output_dir"], "--intermediate-dir", policy["intermediate_dir"],
        "--raw-quality", scope["qc_dir"] + "/raw_read_quality.tsv",
        "--trimmed-quality", policy["intermediate_dir"] + "/terminal_read_quality.tsv",
        "--threads", str(policy["threads"])]


def cohort_result(policy, terminal, manifest):
    directory = ROOT / policy["output_dir"]
    qc = ROOT / terminal["output_dir"]
    result = validate_cohort(directory, manifest, read_table(qc / "run_summary.tsv", []),
                             load_json_yaml(qc / "summary.json"), policy)
    diagnostics = read_table(ROOT / policy["audit_dir"] / "error_diagnostics.tsv", ["read_direction", "last_input_equals_output"])
    if len(diagnostics) != 2 or {row["read_direction"] for row in diagnostics} != {"F", "R"}:
        raise ValueError("Missing error-history audit")
    for row in diagnostics:
        native = result["error_learning"][row["read_direction"]]
        for field in ("self_consistency_rounds", "maximum_rounds"):
            if int(row[field]) != native[field]:
                raise ValueError("Saved model and native error-learning rounds differ")
        for field, native_field in (("repeated_any", "repeated_error_matrix"),
                                    ("converged_before_limit", "converged_before_limit")):
            if row[field] != ("TRUE" if native[native_field] else "FALSE"):
                raise ValueError("Saved model and native convergence flags differ")
    if not all(row["last_input_equals_output"] == "TRUE" for row in diagnostics):
        result["promotion_status"] = "needs_review"
        result["promotion_reason"] = "Final error matrix differs from last input; inspect iteration history"
    result["ecological_interpretation"] = "pending_taxonomy_and_residual_primer_review"
    return result


def audit_products(policy):
    return [*[ROOT / policy["output_dir"] / name for name in DADA_PRODUCTS],
            *[ROOT / policy["audit_dir"] / name for name in AUDIT_PRODUCTS]]


def validate(config_path, policy, terminal, scope, manifest):
    directory = ROOT / policy["audit_dir"]
    expected = [*audit_products(policy), directory / "config_snapshot.yaml", directory / "input_checksums.json"]
    integrity = validate_artifacts(directory, ROOT, expected, policy)
    native = validate_run(ROOT / policy["output_dir"], ROOT)
    cohort = cohort_result(policy, terminal, manifest)
    if cohort != load_json_yaml(directory / "cohort_validation.json"):
        raise ValueError("Cohort validation changed")
    if native != load_json_yaml(ROOT / policy["output_dir"] / "validation.json"):
        raise ValueError("Native DADA2 validation changed")
    return {**integrity, **cohort, "native_validation": native}


def execute(config_path, policy, terminal, scope, manifest, processing):
    audit = ROOT / policy["audit_dir"]
    output = ROOT / policy["output_dir"]
    intermediate = ROOT / policy["intermediate_dir"]
    if (audit / "SUCCESS").exists() or (output / "SUCCESS").exists():
        raise ValueError("Successful execution already exists; use validate instead of overwriting")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Partial execution exists; preserve and review it before a new attempt")
    qc = ROOT / terminal["output_dir"]
    logging.info("Verifying all upstream terminal QC hashes before execution")
    upstream = validate_artifacts(qc, ROOT, [qc / name for name in TERMINAL_PRODUCTS], terminal)
    baseline = load_json_yaml(ROOT / scope["qc_dir"] / "provenance.json")
    raw_quality = ROOT / scope["qc_dir"] / "raw_read_quality.tsv"
    expected_raw = [row for row in baseline["outputs"] if row["path"] == raw_quality.relative_to(ROOT).as_posix()]
    if expected_raw != file_records([raw_quality], ROOT):
        raise ValueError("Raw QC counts differ from the validated baseline")
    terminal_summary = load_json_yaml(qc / "summary.json")
    if terminal_summary["minimum_length_diagnostic"] != processing["dada2"]["filter"]["min_len"]:
        raise ValueError("Terminal length diagnostic differs from the registered filter")
    source_files = inputs(config_path, policy, terminal, scope, manifest)
    fastq_bytes = sum(path.stat().st_size for path in source_files if path.name.endswith(".fastq.gz"))
    if shutil.disk_usage(ROOT).free < policy["minimum_free_space_factor"] * fastq_bytes:
        raise ValueError("Insufficient disk reserve for independent DADA2")
    intermediate.mkdir(parents=True, exist_ok=True)
    quality = quality_rows(manifest, read_table(qc / "run_summary.tsv", []))
    quality_path = intermediate / "terminal_read_quality.tsv"
    write_tsv_atomic(quality_path, quality, list(quality[0]))
    source_files.append(quality_path)
    state = begin_stage(audit, ROOT)
    before = file_records(source_files, ROOT)
    command = dada_command(policy, terminal, scope)
    write_json(intermediate / "execution_start.json", {**state, "inputs": before, "command": command})
    try:
        subprocess.run(["Rscript", "-e",
            'a <- commandArgs(TRUE); stopifnot(as.character(getRversion()) == a[1], as.character(packageVersion("dada2")) == a[2])',
            policy["expected_r_version"], policy["expected_dada2_version"]], cwd=ROOT, check=True)
        subprocess.run([*command, "--validate-only"], cwd=ROOT, check=True)
        logging.info("Executing registered DADA2 with %s threads", policy["threads"])
        with (intermediate / "execution.log").open("w", encoding="utf-8") as log:
            subprocess.run(command, cwd=ROOT, check=True, stdout=log, stderr=subprocess.STDOUT)
        shutil.copyfile(intermediate / "execution.log", audit / "execution_log.txt")
        native = validate_run(output, ROOT)
        write_json(output / "validation.json", native)
        subprocess.run(["Rscript", "scripts/study_dada2/error_diagnostics.R",
            "--intermediate-dir", policy["intermediate_dir"], "--output", str(audit / "error_diagnostics.tsv"),
            "--config", scope["project_config"], "--manifest", scope["manifest"], "--run-dir", policy["output_dir"],
            "--expected-r", policy["expected_r_version"], "--expected-dada2", policy["expected_dada2_version"],
            "--expected-threads", str(policy["threads"])], cwd=ROOT, check=True)
        cohort = cohort_result(policy, terminal, manifest)
        write_json(audit / "cohort_validation.json", cohort)
        finish_stage(audit, ROOT, state, before, config_path, policy, audit_products(policy),
                     upstream, collect_versions("Rscript"), command)
        result = validate(config_path, policy, terminal, scope, manifest)
        write_json(audit / "validation.json", result)
        logging.info("DADA2 integrity validated; promotion status: %s", result["promotion_status"])
    except Exception as error:
        write_json(audit / "failure.json", {**state, "status": "failed", "error": str(error),
                   "execution_log": (intermediate / "execution.log").relative_to(ROOT).as_posix()})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("execute", "validate"))
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    config_path = repository_path(ROOT, str(args.config))
    policy, terminal, scope, manifest, processing = load_policy(config_path)
    if args.action == "execute":
        execute(config_path, policy, terminal, scope, manifest, processing)
    else:
        print(json.dumps(validate(config_path, policy, terminal, scope, manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
