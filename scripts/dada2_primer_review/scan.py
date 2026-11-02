#!/usr/bin/env python3
"""Diagnose full primer-compatible motifs in validated DADA2 ASVs.

Inputs: a root-relative policy referencing DADA2 policy/profile, completed native
and outer validations, ASV FASTA/counts/sequences and their recorded checksums.
Outputs: coordinate.tsv, per_sample.tsv, summary.json and SHA-256 provenance in
the configured separate diagnostic directory. No reads/ASVs are changed or
excluded. Fails on invalid policies, incomplete/stale outputs, inconsistent
counts/FASTA, ambiguous ASV bases or altered diagnostic artifacts. Source FASTQ
are not reopened. Coordinates are zero-based, end-exclusive; matches are not
taxonomic assignments or proof of contamination.
"""
import argparse
from decimal import Decimal, InvalidOperation, ROUND_FLOOR
import importlib.util
import json
import logging
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, read_counts, read_table, unique_rows, validate_artifacts, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.provenance import collect_versions, file_records
from cacao_inventory.taxonomy_reference import repository_path

IUPAC = dict(zip("ACGTRYSWKMBDHVN", ("A", "C", "G", "T", "AG", "CT", "CG", "AT", "GT", "AC", "CGT", "AGT", "ACT", "ACG", "ACGT")))
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
PRODUCTS = ["coordinate.tsv", "per_sample.tsv", "summary.json"]
COORDINATE_COLUMNS = ["study_id", "asv_id", "primer_name", "primer_role", "orientation", "pattern", "start", "end", "errors", "maximum_errors", "sequence_length", "asv_total_reads"]
SAMPLE_COLUMNS = ["study_id", "sample_id", "total_reads", "matched_asvs", "matched_reads", "matched_fraction"]


def maximum_errors(pattern: str, error_rate) -> int:
    if not pattern or any(base not in IUPAC for base in pattern) or not any(base != "N" for base in pattern):
        raise ValueError("Primer must use IUPAC DNA symbols and contain informative bases")
    if isinstance(error_rate, bool):
        raise ValueError("Invalid primer error rate")
    try:
        rate = Decimal(str(error_rate))
    except InvalidOperation as exc:
        raise ValueError("Invalid primer error rate") from exc
    if not rate.is_finite() or not 0 <= rate < 1:
        raise ValueError("Primer error rate must be finite in [0,1)")
    # Cutadapt excludes adapter N from its error-rate denominator; other IUPAC
    # symbols count as one informative position. Read wildcards are disabled.
    return int((rate * sum(base != "N" for base in pattern)).to_integral_value(rounding=ROUND_FLOOR))


def find_matches(sequence: str, pattern: str, error_rate) -> list[dict]:
    if not sequence or not re.fullmatch(r"[ACGT]+", sequence):
        raise ValueError("ASV sequences must be nonempty unambiguous uppercase DNA")
    limit = maximum_errors(pattern, error_rate)
    matches = []
    for start in range(len(sequence) - len(pattern) + 1):
        errors = sum(sequence[start + offset] not in IUPAC[base] for offset, base in enumerate(pattern))
        if errors <= limit:
            matches.append({"start": start, "end": start + len(pattern), "errors": errors})
    return matches


def patterns(processing: dict) -> list[dict]:
    detection = processing["primer_detection"]
    if detection["allow_indels"] is not False or detection["minimum_overlap_fraction"] != 1:
        raise ValueError("Source primer detection must require full matches without indels")
    result = []
    for role in ("forward", "reverse"):
        primer = processing["primers"][role]
        sequence = primer["gene_specific_sequence"]
        maximum_errors(sequence, detection["error_rate"])
        for orientation, pattern in (("as_synthesized", sequence), ("reverse_complement", sequence.translate(COMPLEMENT)[::-1])):
            result.append({"primer_role": role, "primer_name": primer["name"], "orientation": orientation,
                           "pattern": pattern, "maximum_errors": maximum_errors(pattern, detection["error_rate"])})
    return result


def scan_tables(sequences: list[dict], counts: list[dict], processing: dict) -> tuple[list, list, dict]:
    definitions = unique_rows(sequences, ("asv_id",))
    unique_rows(counts, ("study_id", "sample_id"))
    if not definitions or not counts or len({row["study_id"] for row in [*sequences, *counts]}) != 1:
        raise ValueError("ASV tables must describe one nonempty study")
    features = {key[0] for key in definitions}
    totals = {feature: 0 for feature in features}
    for row in counts:
        if set(row) != {"study_id", "sample_id", *features}:
            raise ValueError("ASV count columns differ from sequence definitions")
        for feature in features:
            if not re.fullmatch(r"[0-9]+", str(row[feature])):
                raise ValueError("ASV counts must be nonnegative integers")
            totals[feature] += int(row[feature])
    probes, coordinates, matched = patterns(processing), [], set()
    per_pattern = {(probe["primer_role"], probe["orientation"]): set() for probe in probes}
    for (feature,), row in definitions.items():
        sequence = row["sequence"]
        if str(len(sequence)) != str(row["length"]) or str(totals[feature]) != str(row["total_reads"]) or totals[feature] == 0:
            raise ValueError("ASV length/total differs from sequence or counts")
        for probe in probes:
            for match in find_matches(sequence, probe["pattern"], processing["primer_detection"]["error_rate"]):
                coordinates.append({"study_id": row["study_id"], "asv_id": feature, **probe, **match,
                                    "sequence_length": len(sequence), "asv_total_reads": totals[feature]})
                matched.add(feature)
                per_pattern[probe["primer_role"], probe["orientation"]].add(feature)
    samples = []
    for row in counts:
        total = sum(int(row[feature]) for feature in features)
        if not total:
            raise ValueError("Empty ASV sample is not a completed DADA2 cohort")
        selected = sum(int(row[feature]) for feature in matched)
        samples.append({"study_id": row["study_id"], "sample_id": row["sample_id"], "total_reads": total,
                        "matched_asvs": sum(int(row[feature]) > 0 for feature in matched),
                        "matched_reads": selected, "matched_fraction": selected / total})
    total = sum(totals.values())
    summary = {"study_id": sequences[0]["study_id"], "samples": len(samples), "asvs": len(features),
               "total_reads": total, "matched_asvs": len(matched), "matched_reads": sum(totals[feature] for feature in matched),
               "coordinate_matches": len(coordinates), "asv_exclusions": 0, "sample_exclusions": 0, "sequence_edits": 0,
               "interpretation": "Full primer-compatible motifs only; no taxonomic or contamination identity inferred",
               "patterns": [{**probe, "matched_asvs": len(per_pattern[probe["primer_role"], probe["orientation"]]),
                             "matched_reads": sum(totals[feature] for feature in per_pattern[probe["primer_role"], probe["orientation"]])} for probe in probes]}
    summary["matched_fraction"] = summary["matched_reads"] / total
    return coordinates, samples, summary


def load_policy(config_path: Path, root: Path = ROOT) -> tuple:
    policy = load_json_yaml(config_path)
    expected = {"schema_version": "1.0.0", "method": "full_length_iupac_hamming_all_positions",
                "error_rate_source": "profile.primer_detection.error_rate", "adapter_n_error_denominator": "exclude",
                "read_wildcards": False, "allow_indels": False, "partial_matches": False,
                "trim_sequences": False, "exclude_asvs": False, "sample_exclusions": False}
    if set(policy) != {*expected, "source_dada2_policy", "output_dir"} or any(policy[key] != value or (isinstance(value, bool) and type(policy[key]) is not bool) for key, value in expected.items()):
        raise ValueError("Unsupported primer diagnostic policy")
    dada_path = repository_path(root, policy["source_dada2_policy"])
    dada = load_json_yaml(dada_path)
    terminal_path = repository_path(root, dada["source_terminal_policy"])
    terminal = load_json_yaml(terminal_path)
    scope_path = repository_path(root, terminal["source_scope"])
    scope = load_json_yaml(scope_path)
    profile_path = repository_path(root, scope["project_config"])
    processing = load_json_yaml(profile_path)["amplicon_processing"][scope["bioproject"]]
    output = repository_path(root, policy["output_dir"])
    if not output.is_relative_to((root / "results/diagnostics" / scope["study_id"]).resolve()):
        raise ValueError("Diagnostics must use the study-specific results/diagnostics directory")
    if processing["study_id"] != scope["study_id"]:
        raise ValueError("Profile study differs from source scope")
    patterns(processing)
    return policy, dada, scope, processing, [dada_path, terminal_path, scope_path, profile_path, repository_path(root, scope["manifest"])]


def checked_source(dada: dict, scope: dict, source_paths: list[Path], root: Path = ROOT) -> tuple[dict, list[Path]]:
    directory, audit = (repository_path(root, dada[key]) for key in ("output_dir", "audit_dir"))
    record = load_json_yaml(audit / "provenance.json")
    success = (audit / "SUCCESS").read_text().strip()
    if not re.fullmatch(r"[a-f0-9]{40}", success) or record["status"] != "success" or record["git_commit"] != success or record["parameters"] != dada:
        raise ValueError("Invalid outer DADA2 success/provenance")
    outputs = record["outputs"]
    paths = [repository_path(root, row["path"]) for row in outputs]
    if not paths or len(paths) != len(set(paths)) or any(not (path.is_relative_to(directory) or path.is_relative_to(audit)) or path.suffix.lower() in (".gz", ".fastq", ".fq", ".rds", ".bam") for path in paths):
        raise ValueError("Outer DADA2 outputs must be unique small analysis artifacts")
    if sorted(file_records(paths, root), key=lambda row: row["path"]) != sorted(outputs, key=lambda row: row["path"]):
        raise ValueError("Outer DADA2 output hashes changed")
    bound = {row["path"]: row for row in record["inputs"]}
    if len(bound) != len(record["inputs"]):
        raise ValueError("Duplicate outer DADA2 input records")
    if json.loads((audit / "input_checksums.json").read_text()) != record["inputs"]:
        raise ValueError("Outer DADA2 input checksum list differs from provenance")
    module_path = root / "scripts/dada2/validate_outputs.py"
    for row in file_records([*source_paths, module_path], root):
        if bound.get(row["path"]) != row:
            raise ValueError("Source policy/profile/manifest or native validator changed since DADA2")
    required = [directory / name for name in ("SUCCESS", "validation.json", "asv_counts.tsv", "asv_sequences.tsv", "asv_sequences.fasta", "output_checksums.tsv")]
    required += [audit / "cohort_validation.json", audit / "input_checksums.json"]
    if not set(required).issubset(paths):
        raise ValueError("Outer output hashes do not bind essential DADA2 outputs")
    for row in read_table(directory / "output_checksums.tsv", ["path", "md5"]):
        path = (directory / row["path"]).resolve()
        if not path.is_relative_to(directory) or path not in paths:
            raise ValueError("Native DADA2 output is not bound by outer provenance")
    spec = importlib.util.spec_from_file_location("primer_review_native_validation", module_path)
    native_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(native_module)
    native = native_module.validate_run(directory)  # Deliberately no FASTQ/input rehash.
    saved_native = load_json_yaml(directory / "validation.json")
    if {key: value for key, value in native.items() if key != "inputs_checked"} != {key: value for key, value in saved_native.items() if key != "inputs_checked"}:
        raise ValueError("Native DADA2 validation changed")
    outer = load_json_yaml(audit / "validation.json")
    cohort = load_json_yaml(audit / "cohort_validation.json")
    if outer["status"] != "valid" or outer["source_git_commit"] != success or outer["native_validation"] != saved_native or any(outer.get(key) != value for key, value in cohort.items()):
        raise ValueError("Outer DADA2 validation is incomplete or inconsistent")
    if native["study_id"] != scope["study_id"] or cohort["samples"] != scope["expected_runs"]:
        raise ValueError("DADA2 cohort differs from source scope")
    return {"status": "valid", "source_git_commit": success, "native_validation": native,
            "promotion_status": cohort["promotion_status"], "fastq_rehashed": False}, [*paths, *source_paths, audit / "provenance.json", audit / "SUCCESS", audit / "validation.json", module_path]


def evaluate(dada: dict, processing: dict, root: Path = ROOT) -> tuple:
    directory = repository_path(root, dada["output_dir"])
    sequences = read_table(directory / "asv_sequences.tsv", ["study_id", "asv_id", "sequence", "length", "total_reads"])
    counts, _ = read_counts(directory / "asv_counts.tsv")
    expected_fasta = [line for row in sequences for line in (">" + row["asv_id"], row["sequence"])]
    if (directory / "asv_sequences.fasta").read_text().splitlines() != expected_fasta:
        raise ValueError("ASV FASTA differs from sequence table")
    return scan_tables(sequences, counts, processing)


def run(config_path: Path, action: str):
    policy, dada, scope, processing, source_paths = load_policy(config_path)
    upstream, source_files = checked_source(dada, scope, source_paths)
    output = repository_path(ROOT, policy["output_dir"])
    products = [output / name for name in PRODUCTS]
    expected = [*products, output / "config_snapshot.yaml", output / "input_checksums.json"]
    if action == "validate":
        result = validate_artifacts(output, ROOT, expected, policy)
        coordinates, samples, summary = evaluate(dada, processing)
        for path, rows, columns in ((products[0], coordinates, COORDINATE_COLUMNS), (products[1], samples, SAMPLE_COLUMNS)):
            if read_table(path, columns) != [{key: str(row[key]) for key in columns} for row in rows]:
                raise ValueError("Diagnostic table differs from a fresh deterministic scan")
        if load_json_yaml(output / "summary.json") != summary:
            raise ValueError("Diagnostic summary differs from a fresh deterministic scan")
        return result
    if output.exists() and any(output.iterdir()):
        raise ValueError("Diagnostic directory is not empty; preserve prior results and use validate")
    helpers = [ROOT / "python/cacao_inventory" / (name + ".py") for name in ("analysis_artifacts", "config", "io", "provenance", "taxonomy_reference", "schema")]
    before = file_records([*source_files, config_path, Path(__file__).resolve(), *helpers], ROOT)
    state = begin_stage(output, ROOT)
    coordinates, samples, summary = evaluate(dada, processing)
    write_tsv_atomic(products[0], coordinates, COORDINATE_COLUMNS)
    write_tsv_atomic(products[1], samples, SAMPLE_COLUMNS)
    write_json(products[2], summary)
    finish_stage(output, ROOT, state, before, config_path, policy, products, upstream, collect_versions(), sys.argv)
    checked = validate_artifacts(output, ROOT, expected, policy)
    write_json(output / "validation.json", {**checked, **summary})
    return {**checked, **summary}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("execute", "validate"))
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    logging.info("Primer review: %s", run(repository_path(ROOT, str(args.config)), args.action))


if __name__ == "__main__":
    main()
