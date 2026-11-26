#!/usr/bin/env python3
"""Execute or validate one registered study taxonomy without exclusions.

Inputs: execution policy, immutable taxonomy profile/reference, validated DADA2
and a validated full-ASV primer review. Outputs: native taxonomy tables/figures
and a separate SHA-256 audit with execution log. Fails on changed identities,
counts, metadata, source hashes, software/settings, incomplete primer review,
unready DADA2, overlapping paths, or any complete/partial existing execution.
Only the native classifier performs classification; the wrapper does not retry,
change thresholds or discard flagged ASVs. Validation never classifies again.
"""
import argparse
from datetime import datetime
import importlib.util
import json
import logging
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/taxonomy"))
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, read_table, validate_artifacts, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.provenance import collect_versions, file_records, sha256_file
from cacao_inventory.taxonomy_reference import repository_path, validate_taxonomy_config
from seeded_policy import validate_native_rng


def module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


REVIEW = module(ROOT / "scripts/dada2_primer_review/scan.py", "study_taxonomy_primer_review")
NATIVE = module(ROOT / "scripts/taxonomy/validate_taxonomy.py", "study_taxonomy_native_validator")
REPETITION = module(ROOT / "scripts/taxonomy/audit_repetition.py", "study_taxonomy_repetition")
TAXONOMY_PRODUCTS = sorted(NATIVE.REQUIRED | {"SUCCESS", "provenance.json", "validation.json"})
AUDIT_PRODUCTS = ["source_validation.json", "execution_log.txt", "repeat_execution_log.txt", "repetition.json"]
AUDIT_METADATA = ["SUCCESS", "validation.json", "provenance.json", "config_snapshot.yaml", "input_checksums.json"]
HELPERS = ("analysis_artifacts", "config", "download", "io", "provenance", "schema", "taxonomy_reference", "taxonomy_validation")


def load_policy(config_path: Path, root: Path = ROOT) -> tuple:
    policy = load_json_yaml(config_path)
    fixed = {"schema_version": "1.0.0", "method": "independent_study_registered_taxonomy",
             "require_dada2_promotion": "ready_for_taxonomy", "residual_primer_policy": "retain_all_diagnosed_asvs",
             "sample_exclusions": False, "asv_exclusions": False, "parameter_overrides": False, "processing_ceiling": "taxonomy",
             "require_identical_repetition": True}
    variable = {"source_dada2_policy", "source_primer_review_policy", "taxonomy_profile", "taxonomy_profile_sha256",
                "output_dir", "audit_dir", "repeat_dir", "repeat_threads", "expected_samples", "expected_asvs", "expected_nonchim_pairs", "threads", "memory_mb"}
    if set(policy) != set(fixed) | variable or any(type(policy[key]) is not type(value) or policy[key] != value for key, value in fixed.items()):
        raise ValueError("Unsupported independent taxonomy policy")
    review_path = repository_path(root, policy["source_primer_review_policy"])
    review, dada, scope, processing, source_paths = REVIEW.load_policy(review_path, root)
    if review["source_dada2_policy"] != policy["source_dada2_policy"]:
        raise ValueError("Taxonomy and primer review reference different DADA2 policies")
    profile_path = repository_path(root, policy["taxonomy_profile"])
    if not re.fullmatch(r"[a-f0-9]{64}", policy["taxonomy_profile_sha256"]) or sha256_file(profile_path) != policy["taxonomy_profile_sha256"]:
        raise ValueError("Registered taxonomy profile hash changed")
    profile = load_json_yaml(profile_path)
    validate_taxonomy_config(profile, root)
    validate_native_rng(profile)
    classifier = profile["classification"]
    if (not classifier.get("native_rng") or classifier["native_rng"]["dada2_version"] != dada["expected_dada2_version"]
            or classifier["random_seed"] != processing["dada2"]["random_seed"]):
        raise ValueError("Taxonomy adapter/version/seed differs from the registered study")
    for field in ("expected_samples", "expected_asvs", "expected_nonchim_pairs", "threads", "repeat_threads", "memory_mb"):
        if type(policy[field]) is not int or policy[field] <= 0:
            raise ValueError(f"Invalid positive integer: {field}")
    if (policy["expected_samples"] != dada["expected_runs"] or policy["expected_samples"] != scope["expected_runs"]
            or policy["threads"] != classifier["threads"] or policy["repeat_threads"] >= policy["threads"]
            or policy["memory_mb"] > scope["resources"]["memory_budget_mb"]):
        raise ValueError("Taxonomy cohort/thread/memory policy disagrees with registered settings")
    output, audit = [repository_path(root, policy[key]) for key in ("output_dir", "audit_dir")]
    repeat = repository_path(root, policy["repeat_dir"])
    parent = (root / "results/taxonomy" / scope["study_id"]).resolve()
    protected = [repository_path(root, scope["taxonomy_dir"]), repository_path(root, dada["output_dir"]),
                 repository_path(root, dada["audit_dir"]), repository_path(root, review["output_dir"])]
    if any(path == parent or not path.is_relative_to(parent) for path in (output, audit)) or output.parent != audit.parent:
        raise ValueError("Taxonomy output and audit must be sibling study directories")
    repeat_parent = (root / "results/intermediate/taxonomy" / scope["study_id"]).resolve()
    if repeat == repeat_parent or not repeat.is_relative_to(repeat_parent):
        raise ValueError("Taxonomy repetition must use the study's ignored intermediate directory")
    for index, path in enumerate((output, audit, repeat)):
        if any(path == other or path in other.parents or other in path.parents for other in [*(output, audit, repeat)[index + 1:], *protected]):
            raise ValueError("Taxonomy paths overlap source, baseline or one another")
    return policy, profile, review, dada, scope, processing, source_paths


def inputs(config_path: Path, root: Path = ROOT) -> list[Path]:
    """List reproducibility inputs for the workflow without classification."""
    policy, profile, review, dada, scope, _, source_paths = load_policy(config_path, root)
    source_audit = repository_path(root, dada["audit_dir"])
    source_records = load_json_yaml(source_audit / "provenance.json")["outputs"]
    primer_dir = repository_path(root, review["output_dir"])
    terminal = load_json_yaml(repository_path(root, dada["source_terminal_policy"]))
    result = [Path(config_path).resolve(), *source_paths, repository_path(root, policy["taxonomy_profile"]),
              repository_path(root, policy["source_primer_review_policy"]), root / "config/config.yaml",
              *[repository_path(root, row["path"]) for row in source_records],
              *[source_audit / name for name in ("SUCCESS", "validation.json", "provenance.json")],
              *[primer_dir / name for name in (*REVIEW.PRODUCTS, "SUCCESS", "validation.json", "provenance.json", "config_snapshot.yaml", "input_checksums.json")],
              *[repository_path(root, terminal["output_dir"]) / name for name in ("run_summary.tsv", "summary.json")],
              repository_path(root, profile["reference"]["path"]), repository_path(root, profile["reference"]["provenance_path"]),
              root / "scripts/study_taxonomy/execute.py", root / "scripts/dada2_primer_review/scan.py",
              root / "scripts/study_dada2/cohort.py", root / "scripts/dada2/helpers.R", root / "scripts/dada2/validate_outputs.py",
              root / "environment/conda-linux-64.lock", root / "workflow/study_taxonomy.Snakefile",
              root / "workflow/Snakefile", root / "workflow/study.Snakefile", root / "workflow/rules/pilot_taxonomy.smk",
              *[root / "python/cacao_inventory" / f"{name}.py" for name in HELPERS],
              *[path for path in (root / "scripts/taxonomy").rglob("*") if path.suffix in (".py", ".R", ".cpp")]]
    return sorted(set(path.resolve() for path in result))


def check_sources(context: tuple, root: Path = ROOT) -> dict:
    policy, _, review, dada, scope, processing, source_paths = context
    upstream, _ = REVIEW.checked_source(dada, scope, source_paths, root)
    if upstream["promotion_status"] != policy["require_dada2_promotion"]:
        raise ValueError("DADA2 has not passed the registered taxonomy promotion gate")
    native = upstream["native_validation"]
    for field, expected in (("samples", "expected_samples"), ("asvs", "expected_asvs"), ("nonchim_pairs", "expected_nonchim_pairs")):
        if native[field] != policy[expected]:
            raise ValueError(f"DADA2 source differs from registered {expected}")
    source_audit = repository_path(root, dada["audit_dir"])
    terminal = load_json_yaml(repository_path(root, dada["source_terminal_policy"]))
    terminal_qc = repository_path(root, terminal["output_dir"])
    cohort_path = root / "scripts/study_dada2/cohort.py"
    bound = {row["path"]: row for row in load_json_yaml(source_audit / "provenance.json")["inputs"]}
    for row in file_records([cohort_path, terminal_qc / "run_summary.tsv", terminal_qc / "summary.json"], root):
        if bound.get(row["path"]) != row:
            raise ValueError("Cohort validator or preceding QC changed since DADA2")
    cohort_module = module(cohort_path, "study_taxonomy_source_cohort")
    cohort = cohort_module.validate_cohort(repository_path(root, dada["output_dir"]), read_table(repository_path(root, scope["manifest"]), []),
                                         read_table(terminal_qc / "run_summary.tsv", []), load_json_yaml(terminal_qc / "summary.json"), dada)
    saved_cohort = load_json_yaml(source_audit / "cohort_validation.json")
    if any(saved_cohort.get(key) != value for key, value in cohort.items()):
        raise ValueError("DADA2 cohort metadata/counts or promotion state changed")
    primer_dir = repository_path(root, review["output_dir"])
    integrity = validate_artifacts(primer_dir, root, [primer_dir / name for name in (*REVIEW.PRODUCTS, "config_snapshot.yaml", "input_checksums.json")], review)
    coordinates, samples, summary = REVIEW.evaluate(dada, processing, root)
    for name, rows, columns in (("coordinate.tsv", coordinates, REVIEW.COORDINATE_COLUMNS), ("per_sample.tsv", samples, REVIEW.SAMPLE_COLUMNS)):
        if read_table(primer_dir / name, columns) != [{key: str(row[key]) for key in columns} for row in rows]:
            raise ValueError("Primer diagnostic table disagrees with unchanged source ASVs")
    if load_json_yaml(primer_dir / "summary.json") != summary or load_json_yaml(primer_dir / "validation.json") != {**integrity, **summary}:
        raise ValueError("Primer diagnostic summary/validation is inconsistent")
    return {"status": "valid", "dada2": upstream, "cohort": cohort, "primer_review": {**integrity, **summary}}


def native_result(context: tuple, root: Path = ROOT, repeat: bool = False, expected_revision: str | None = None) -> dict:
    policy, profile, _, dada, scope, _, _ = context
    output = repository_path(root, policy["repeat_dir"] if repeat else policy["output_dir"])
    threads = policy["repeat_threads"] if repeat else policy["threads"]
    result = NATIVE.validate(output, repository_path(root, dada["output_dir"]), root)
    expected = {"study_id": scope["study_id"], "samples": policy["expected_samples"], "asvs": policy["expected_asvs"],
                "total_reads": policy["expected_nonchim_pairs"], "excluded_asvs": 0, "status": "valid"}
    if any(result.get(key) != value for key, value in expected.items()):
        raise ValueError("Native taxonomy validation changed cohort, counts or exclusions")
    if load_json_yaml(output / "config_snapshot.yaml") != profile:
        raise ValueError("Executed taxonomy profile differs from the immutable profile")
    provenance = load_json_yaml(output / "provenance.json")
    if provenance["effective_threads"] != threads:
        raise ValueError("Executed taxonomy thread count differs from policy")
    if expected_revision is not None and provenance["git_commit"] != expected_revision:
        raise ValueError("Native taxonomy producer differs from the outer execution")
    if (provenance["classification"] != profile["classification"] or provenance["reference"] != profile["reference"]
            or provenance["rng_kind"] != profile["classification"]["rng_kind"]):
        raise ValueError("Native taxonomy settings or seed differ from the immutable profile")
    expected_command = ["Rscript", "scripts/taxonomy/assign_taxonomy.R", "--config", str(repository_path(root, policy["taxonomy_profile"])),
                        "--input-dir", str(repository_path(root, dada["output_dir"])), "--output-dir", str(output), "--threads", str(threads)]
    if provenance["command"] != expected_command:
        raise ValueError("Native taxonomy command differs from requested source/profile/output")
    input_paths = {repository_path(root, row["path"]) for row in provenance["inputs"]}
    dada_directory = repository_path(root, dada["output_dir"])
    expected_dada = {dada_directory / name for name in ("SUCCESS", "validation.json", "asv_sequences.tsv", "asv_counts.tsv", "sample_metadata.tsv")}
    observed_dada = {path for path in input_paths if path.is_relative_to((root / "results/dada2").resolve())}
    required = {repository_path(root, policy["taxonomy_profile"]), root / "config/config.yaml",
                repository_path(root, profile["reference"]["path"]), repository_path(root, profile["reference"]["provenance_path"]), *expected_dada}
    if observed_dada != expected_dada or not required.issubset(input_paths):
        raise ValueError("Native taxonomy inputs do not bind the requested DADA2 cohort/profile/reference")
    if provenance["upstream_validation"] != load_json_yaml(dada_directory / "validation.json"):
        raise ValueError("Native taxonomy upstream validation differs from the requested DADA2 execution")
    versions = read_table(output / "software_versions.tsv", ["software", "version"])
    observed = {row["software"]: row["version"] for row in versions}
    expected_versions = {"R": dada["expected_r_version"], "dada2": dada["expected_dada2_version"],
                         "taxonomy_native_rng": profile["classification"]["native_rng"]["adapter"]}
    if len(observed) != len(versions) or any(observed.get(key) != value for key, value in expected_versions.items()):
        raise ValueError("Executed taxonomy software/seeded adapter differs from registration")
    return result


def products(policy: dict, root: Path = ROOT) -> list[Path]:
    return [*[repository_path(root, policy[key]) / name for key in ("output_dir", "repeat_dir") for name in TAXONOMY_PRODUCTS],
            *[repository_path(root, policy["audit_dir"]) / name for name in AUDIT_PRODUCTS]]


def check_repetition(context: tuple, root: Path = ROOT) -> dict:
    policy, profile = context[:2]
    baseline, repeat = (repository_path(root, policy[key]) for key in ("output_dir", "repeat_dir"))
    record = load_json_yaml(repository_path(root, policy["audit_dir"]) / "repetition.json")
    tables = REPETITION.compare_tables(baseline, repeat)
    if record["all_identical"] is not True or record["tables"] != tables or not all(row["identical"] for row in tables.values()):
        raise ValueError("Taxonomy repetition differs from the primary classification or saved comparison")
    histories = {}
    for label, directory, threads in (("baseline", baseline, policy["threads"]), ("repeat", repeat, policy["repeat_threads"])):
        provenance = load_json_yaml(directory / "provenance.json")
        saved = record["executions"][label]
        if ((root / saved["directory"]).resolve() != directory or saved["config"] != profile
                or saved["git_commit"] != (directory / "SUCCESS").read_text().strip() or saved["effective_threads"] != threads
                or saved["provenance_sha256"] != sha256_file(directory / "provenance.json")
                or any(saved[key] != provenance[key] for key in ("started_at_utc", "finished_at_utc"))):
            raise ValueError("Taxonomy repetition execution metadata differ from recorded runs")
        histories[label] = provenance
    if sorted(histories["baseline"]["inputs"], key=lambda row: row["path"]) != sorted(histories["repeat"]["inputs"], key=lambda row: row["path"]):
        raise ValueError("Taxonomy repetitions used different input bytes")
    begin = datetime.fromisoformat(histories["baseline"]["started_at_utc"])
    finish = datetime.fromisoformat(histories["baseline"]["finished_at_utc"])
    repeat_begin = datetime.fromisoformat(histories["repeat"]["started_at_utc"])
    repeat_finish = datetime.fromisoformat(histories["repeat"]["finished_at_utc"])
    if any(value.tzinfo is None for value in (begin, finish, repeat_begin, repeat_finish)) or not begin <= finish <= repeat_begin <= repeat_finish:
        raise ValueError("Taxonomy primary/repetition chronology is inconsistent")
    return {"all_identical": True, "tables_checked": len(tables), "primary_threads": policy["threads"], "repeat_threads": policy["repeat_threads"]}


def validate(config_path: Path, root: Path = ROOT) -> dict:
    context = load_policy(config_path, root)
    policy = context[0]
    audit = repository_path(root, policy["audit_dir"])
    integrity = validate_artifacts(audit, root, [*products(policy, root), audit / "config_snapshot.yaml", audit / "input_checksums.json"], policy)
    recorded_inputs = load_json_yaml(audit / "provenance.json")["inputs"]
    if {row["path"] for row in recorded_inputs} != {path.relative_to(root).as_posix() for path in inputs(config_path, root)}:
        raise ValueError("Outer taxonomy provenance omits or substitutes required inputs")
    source = check_sources(context, root)
    native = native_result(context, root, expected_revision=integrity["source_git_commit"])
    repeated = native_result(context, root, repeat=True, expected_revision=integrity["source_git_commit"])
    if load_json_yaml(audit / "source_validation.json") != source or load_json_yaml(repository_path(root, policy["output_dir"]) / "validation.json") != native:
        raise ValueError("Stored taxonomy validation no longer matches independent validation")
    if load_json_yaml(repository_path(root, policy["repeat_dir"]) / "validation.json") != repeated:
        raise ValueError("Stored repetition validation no longer matches independent validation")
    repetition = check_repetition(context, root)
    return {**integrity, **{key: value for key, value in native.items() if key not in ("source_git_commit", "artifacts_checked", "inputs_checked")},
            "native_validation": native, "primer_matched_asvs_retained": source["primer_review"]["matched_asvs"],
            "primer_matched_reads_retained": source["primer_review"]["matched_reads"], "sample_exclusions": 0, "repetition": repetition}


def execute(config_path: Path, root: Path = ROOT) -> dict:
    context = load_policy(config_path, root)
    policy, _, _, dada, _, _, _ = context
    output, audit, repeat = (repository_path(root, policy[key]) for key in ("output_dir", "audit_dir", "repeat_dir"))
    if any(path.exists() and (not path.is_dir() or any(path.iterdir())) for path in (output, audit, repeat)):
        raise ValueError("Existing complete or partial taxonomy execution must be preserved; use validate")
    source = check_sources(context, root)
    before = file_records(inputs(config_path, root), root)
    state = begin_stage(audit, root)
    command = [sys.executable, "scripts/taxonomy/run_taxonomy.py", "--config", policy["taxonomy_profile"],
               "--input-dir", dada["output_dir"], "--output-dir", policy["output_dir"], "--threads", str(policy["threads"])]
    try:
        with (audit / "execution_log.txt").open("w", encoding="utf-8") as log:
            subprocess.run(["Rscript", "-e", 'a <- commandArgs(TRUE); stopifnot(as.character(getRversion()) == a[1], as.character(packageVersion("dada2")) == a[2])',
                            dada["expected_r_version"], dada["expected_dada2_version"]], cwd=root, check=True, stdout=log, stderr=subprocess.STDOUT)
            subprocess.run(command, cwd=root, check=True, stdout=log, stderr=subprocess.STDOUT)
        native = native_result(context, root, expected_revision=state["git_commit"])
        write_json(output / "validation.json", native)
        repeat_command = [sys.executable, "scripts/taxonomy/run_taxonomy.py", "--config", policy["taxonomy_profile"],
                          "--input-dir", dada["output_dir"], "--output-dir", policy["repeat_dir"], "--threads", str(policy["repeat_threads"])]
        with (audit / "repeat_execution_log.txt").open("w", encoding="utf-8") as log:
            subprocess.run(repeat_command, cwd=root, check=True, stdout=log, stderr=subprocess.STDOUT)
        repeated = native_result(context, root, repeat=True, expected_revision=state["git_commit"])
        write_json(repeat / "validation.json", repeated)
        subprocess.run([sys.executable, "scripts/taxonomy/audit_repetition.py", "--baseline", policy["output_dir"],
                        "--repeat", policy["repeat_dir"], "--output", str(audit / "repetition.json"), "--require-identical"],
                       cwd=root, check=True, capture_output=True, text=True)
        check_repetition(context, root)
        write_json(audit / "source_validation.json", source)
        finish_stage(audit, root, state, before, config_path, policy, products(policy, root), source, collect_versions("Rscript"),
                     [*command, "THEN", *repeat_command])
        result = validate(config_path, root)
        write_json(audit / "validation.json", result)
        return result
    except Exception as error:
        (audit / "SUCCESS").unlink(missing_ok=True)
        (audit / "validation.json").unlink(missing_ok=True)
        write_json(audit / "failure.json", {**state, "status": "failed", "error": str(error)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("execute", "validate"))
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    path = repository_path(ROOT, str(args.config))
    result = execute(path) if args.action == "execute" else validate(path)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
