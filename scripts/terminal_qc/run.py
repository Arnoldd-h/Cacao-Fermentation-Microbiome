#!/usr/bin/env python3
"""Conservative terminal-primer cleanup and full-cohort QC.

Inputs: registered policy, validated baseline QC/manifest and 5-prime FASTQ.
Outputs: separate paired FASTQ, per-run audits, aggregate QC and SHA-256 record.
Fails on changed baseline hashes, unsupported methods, altered read identities,
nonterminal edits, count inconsistencies or incomplete quality reports.
No read/sample filtering; the source DADA2 minimum length is diagnostic only.
"""

import argparse
from collections import Counter
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, read_table, validate_artifacts, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.pilot_workflow import load_pilot_manifest, processing_configuration
from cacao_inventory.provenance import collect_versions, file_records
from cacao_inventory.qc import cutadapt_detection_counts, reverse_complement_iupac
from cacao_inventory.study_qc import trimmed_path
from cacao_inventory.study_scope import validate_scope
from cacao_inventory.taxonomy_reference import repository_path
from audit import audit_pair, quality_coverage, read_fastqc

PRODUCTS = ["run_summary.tsv", "read_length_counts.tsv", "quality_coverage.tsv", "summary.json"]
REPORTS = ["cutadapt.json", "R1_terminal.json", "R1_any_position.json", "R2_terminal.json", "R2_any_position.json"]
HELPERS = ("analysis_artifacts", "config", "download", "io", "pilot_workflow", "provenance", "qc", "schema", "study_qc", "study_scope", "taxonomy_reference")


def load_policy(path, root=ROOT):
    policy = load_json_yaml(path)
    expected = {"schema_version": "1.0.0", "method": "opposite_primer_rc_full_terminal_single_pass",
                "allow_indels": False, "quality_trimming": False, "discard_reads": False,
                "partial_matches": False, "internal_matches": False, "passes": 1,
                "processing_ceiling": "qc", "length_diagnostic": "source_dada2_min_len_without_filtering"}
    if any(type(policy.get(k)) is not type(v) or policy[k] != v for k, v in expected.items()):
        raise ValueError("Unsupported terminal cleanup method")
    scope = load_json_yaml(repository_path(root, policy["source_scope"]))
    validate_scope(scope, root)
    project = load_json_yaml(root / scope["project_config"])
    manifest = load_pilot_manifest(root / scope["manifest"])
    if len(manifest) != scope["expected_runs"] or any(row["study_id"] != scope["study_id"] or row["bioproject"] != scope["bioproject"] for row in manifest):
        raise ValueError("Terminal cleanup requires the complete source manifest")
    processing = processing_configuration(project, manifest[0])
    if policy["error_rate"] != processing["primer_detection"]["error_rate"]:
        raise ValueError("Terminal cleanup must retain the registered screening tolerance")
    paths = {k: repository_path(root, policy[k]) for k in ("interim_dir", "intermediate_dir", "output_dir")}
    for key, parent in (("interim_dir", "data/interim"), ("intermediate_dir", "results/intermediate/qc"), ("output_dir", "results/qc")):
        study_root = root / parent / scope["study_id"]
        if not paths[key].is_relative_to(study_root) or paths[key] == study_root:
            raise ValueError("Invalid terminal output namespace")
    baseline = [(root / scope[k]).resolve() for k in ("interim_dir", "qc_dir")]
    if any(a == b or a in b.parents or b in a.parents for a in paths.values() for b in baseline):
        raise ValueError("Terminal outputs must not overlap baseline directories")
    if len(set(paths.values())) != 3 or policy["minimum_free_space_factor"] < 2:
        raise ValueError("Invalid output separation or disk reserve")
    return policy, scope, manifest, processing


def method_inputs(config_path, policy, scope):
    return [config_path, ROOT / policy["source_scope"], ROOT / scope["project_config"], ROOT / scope["manifest"],
            ROOT / "environment/conda-linux-64.lock", ROOT / "workflow/terminal_qc.Snakefile",
            *list((ROOT / "scripts/terminal_qc").glob("*.py")),
            *[ROOT / "python/cacao_inventory" / (name + ".py") for name in HELPERS]]


def sequences(processing):
    return [reverse_complement_iupac(processing["primers"][role]["gene_specific_sequence"])
            for role in ("reverse", "forward")]


def read_paths(policy, scope, row):
    inputs = [ROOT / trimmed_path(scope, row, str(d)) for d in (1, 2)]
    outputs = [ROOT / policy["interim_dir"] / row["run_accession"] / f'{row["run_accession"]}_{d}.fastq.gz' for d in (1, 2)]
    return inputs, outputs


def cutadapt_command(policy, primers, inputs, outputs, report):
    # Anchoring requires a full suffix; Cutadapt rejects min_overlap on anchors.
    return ["cutadapt", "--cores", "1", "--no-indels", "--error-rate", str(policy["error_rate"]),
            "--times", "1", "-a", "opposite_rc=" + primers[0] + "$", "-A", "opposite_rc=" + primers[1] + "$",
            "--json", str(report), "-o", str(outputs[0]), "-p", str(outputs[1]), *map(str, inputs)]


def source_validation(scope):
    directory = ROOT / scope["qc_dir"]
    validation = load_json_yaml(directory / "validation.json")
    if validation["status"] != "valid" or validation["study_id"] != scope["study_id"]:
        raise ValueError("A validated baseline QC for this study is required")
    provenance = load_json_yaml(directory / "provenance.json")
    records = {r["path"]: r for r in provenance["inputs"]}
    for record in file_records([ROOT / scope["project_config"], ROOT / scope["manifest"]], ROOT):
        if records.get(record["path"]) != record:
            raise ValueError("Source processing profile or manifest changed since baseline validation")
    return validation, records


def screen_read(path, primer, policy, report, terminal):
    command = ["cutadapt", "--action=none", "--no-indels", "--error-rate", str(policy["error_rate"])]
    command += ["-a", "probe=" + primer + "$"] if terminal else ["--overlap", str(len(primer)), "--front", "probe=" + primer]
    command += ["--json", str(report), "-o", os.devnull, str(path)]
    subprocess.run(command, cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    result = load_json_yaml(report)
    if result["cutadapt_version"] != policy["cutadapt_version"]:
        raise ValueError("Cutadapt version changed")
    return cutadapt_detection_counts(result)


def validate_run_record(audit, directory, method, policy):
    """Reject stale methods and missing/modified per-run execution evidence."""
    if audit["method_checksums"] != method or audit["cutadapt_version"] != policy["cutadapt_version"]:
        raise ValueError("Per-run audit was produced by a different method")
    reports = [directory / name for name in REPORTS]
    if file_records(reports, ROOT) != audit["report_checksums"]:
        raise ValueError("Required per-run Cutadapt reports changed or disappeared")
    cut = load_json_yaml(directory / "cutadapt.json")
    if cut["read_counts"]["input"] != audit["pairs"] or cut["read_counts"]["output"] != audit["pairs"]:
        raise ValueError("Run report and audit disagree on pairs")
    for i, direction in enumerate(("R1", "R2"), 1):
        if cut["read_counts"][f"read{i}_with_adapter"] != audit["directions"][direction]["changed_reads"]:
            raise ValueError("Run report and audit disagree on trimmed reads")
        for label in ("terminal", "any_position"):
            result = load_json_yaml(directory / f"{direction}_{label}.json")
            if cutadapt_detection_counts(result) != (audit["pairs"], audit["directions"][direction]["residual_" + label]):
                raise ValueError("Run screen and audit disagree")
    return reports


def trim_run(config_path, policy, scope, processing, row):
    inputs, outputs = read_paths(policy, scope, row)
    _, source_records = source_validation(scope)
    before = file_records(inputs, ROOT)
    method = file_records(method_inputs(config_path, policy, scope), ROOT)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    if any(source_records.get(record["path"]) != record for record in before):
        raise ValueError("Baseline FASTQ no longer match their validated hashes")
    if shutil.disk_usage(ROOT).free < sum(p.stat().st_size for p in inputs) * policy["minimum_free_space_factor"]:
        raise ValueError("Insufficient free space for terminal cleanup")
    directory = ROOT / policy["intermediate_dir"] / row["run_accession"]
    directory.mkdir(parents=True, exist_ok=True)
    outputs[0].parent.mkdir(parents=True, exist_ok=True)
    temporary = [p.with_name(p.name.replace(".fastq.gz", ".partial.fastq.gz")) for p in outputs]
    primers = sequences(processing)
    report = directory / "cutadapt.json"
    command = cutadapt_command(policy, primers, inputs, temporary, report)
    subprocess.run(command, cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    cut = load_json_yaml(report)
    audit = audit_pair(inputs, temporary, list(map(len, primers)), processing["dada2"]["filter"]["min_len"])
    if cut["cutadapt_version"] != policy["cutadapt_version"] or cut["command_line_arguments"] != command[1:]:
        raise ValueError("Cutadapt version/command differs from the registered operation")
    if cut["read_counts"]["input"] != audit["pairs"] or cut["read_counts"]["output"] != audit["pairs"]:
        raise ValueError("Cutadapt and independent pair counts disagree")
    if any(value is not None for value in cut["read_counts"]["filtered"].values()) or cut["basepair_counts"]["quality_trimmed"] is not None:
        raise ValueError("Unregistered filtering or quality trimming")
    for i, direction in enumerate(("R1", "R2"), 1):
        measured = audit["directions"][direction]
        if measured["changed_reads"] != cut["read_counts"][f"read{i}_with_adapter"]:
            raise ValueError("Cutadapt and independently changed reads disagree")
        for stage, field in (("input", "baseline_length_counts"), ("output", "length_counts")):
            bases = sum(int(length) * count for length, count in measured[field].items())
            if bases != cut["basepair_counts"][f"{stage}_read{i}"]:
                raise ValueError("Cutadapt and independently counted bases disagree")
        for label, terminal in (("terminal", True), ("any_position", False)):
            examined, matched = screen_read(temporary[i-1], primers[i-1], policy, directory / f"{direction}_{label}.json", terminal)
            if examined != audit["pairs"]:
                raise ValueError("Residual screen and pair counts disagree")
            measured["residual_" + label] = matched
    if file_records(inputs, ROOT) != before:
        raise ValueError("Baseline changed during cleanup")
    if file_records(method_inputs(config_path, policy, scope), ROOT) != method:
        raise ValueError("Method changed during cleanup")
    for source, destination in zip(temporary, outputs):
        source.replace(destination)
    write_json(directory / "audit.json", {**audit, "study_id": row["study_id"], "run_accession": row["run_accession"],
               "sample_id": row["sample_id"], "input_checksums": before, "output_checksums": file_records(outputs, ROOT),
               "method_checksums": method, "producer_git_commit": revision,
               "report_checksums": file_records([directory / name for name in REPORTS], ROOT),
               "cutadapt_version": cut["cutadapt_version"]})
    logging.info("%s: %d pairs conserved; terminal removals R1=%d, R2=%d", row["run_accession"], audit["pairs"],
                 audit["directions"]["R1"]["changed_reads"], audit["directions"]["R2"]["changed_reads"])


def finalize(config_path, policy, scope, manifest, processing):
    output = ROOT / policy["output_dir"]
    state = begin_stage(output, ROOT)
    validation, _ = source_validation(scope)
    inputs = [*method_inputs(config_path, policy, scope),
              ROOT / scope["qc_dir"] / "validation.json", ROOT / scope["qc_dir"] / "provenance.json",
              ]
    method = file_records(method_inputs(config_path, policy, scope), ROOT)
    producer_revisions = set()
    summaries, lengths, qualities = [], [], []
    totals = Counter()
    for row in manifest:
        run = row["run_accession"]
        directory = ROOT / policy["intermediate_dir"] / run
        audit = load_json_yaml(directory / "audit.json")
        if (audit["study_id"], audit["sample_id"], audit["run_accession"]) != (row["study_id"], row["sample_id"], run):
            raise ValueError("Audit identity mismatch")
        reports = validate_run_record(audit, directory, method, policy)
        producer_revisions.add(audit["producer_git_commit"])
        reads, cleaned = read_paths(policy, scope, row)
        inputs += [*reads, *cleaned, directory / "audit.json", *reports]
        if file_records(reads, ROOT) != audit["input_checksums"] or file_records(cleaned, ROOT) != audit["output_checksums"]:
            raise ValueError("Audited reads changed")
        for field in ("pairs", "pairs_below_min_len", "pairs_below_min_len_baseline", "pairs_newly_below_min_len"):
            totals[field] += audit[field]
        for i, direction in enumerate(("R1", "R2"), 1):
            item = audit["directions"][direction]
            base = {"study_id": row["study_id"], "sample_id": row["sample_id"], "run_accession": run, "direction": direction}
            summary = {**base, "pairs": audit["pairs"], **{k: v for k, v in item.items() if not k.endswith("length_counts")}}
            for stage, hist in (("baseline", item["baseline_length_counts"]), ("terminal", item["length_counts"])):
                lengths += [{**base, "stage": stage, "length": length, "reads": count} for length, count in sorted(hist.items(), key=lambda kv: int(kv[0]))]
                fastqc = ((ROOT / scope["qc_dir"] / "fastqc_trimmed") if stage == "baseline" else
                          ROOT / policy["intermediate_dir"] / "fastqc") / f"{run}_{i}_fastqc.zip"
                inputs.append(fastqc)
                report = read_fastqc(fastqc)
                if int(report["basic"]["Total Sequences"]) != audit["pairs"] or report["basic"]["Filename"] != cleaned[i-1].name:
                    raise ValueError("FastQC identity/count mismatch")
                status = report["module_status"]["Per base sequence quality"]
                summary[stage + "_base_quality_status"] = status
                totals[stage + "_" + direction + "_base_quality_" + status] += 1
                qualities += [{**base, "stage": stage, **record} for record in quality_coverage(report["quality_rows"], hist)]
            summaries.append(summary)
            for key in ("changed_reads", "bases_removed", "residual_terminal", "residual_any_position", "newly_below_min_len"):
                totals[direction + "_" + key] += item[key]
    if totals["pairs"] != validation["trimmed_pairs"]:
        raise ValueError("Cohort pair counts differ from baseline")
    before = file_records(inputs, ROOT)
    summary = {"study_id": scope["study_id"], "runs": len(manifest), "sample_exclusions": 0,
               "paired_retention_fraction": 1.0, "minimum_length_diagnostic": processing["dada2"]["filter"]["min_len"],
               "length_filter_applied": False, "processing_ceiling": "qc", "producer_git_commits": sorted(producer_revisions), **dict(totals)}
    for name, rows in (("run_summary.tsv", summaries), ("read_length_counts.tsv", lengths), ("quality_coverage.tsv", qualities)):
        write_tsv_atomic(output / name, rows, list(rows[0]))
    write_json(output / "summary.json", summary)
    products = [output / name for name in PRODUCTS]
    finish_stage(output, ROOT, state, before, config_path, policy, products, validation, collect_versions("cutadapt", "fastqc"), sys.argv)
    checked = validate_artifacts(output, ROOT, [*products, output / "config_snapshot.yaml", output / "input_checksums.json"], policy)
    write_json(output / "validation.json", {**checked, **summary})
    logging.info("Terminal QC complete: %s", summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("trim", "finalize", "validate"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-accession")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    policy, scope, manifest, processing = load_policy(args.config)
    if args.step == "trim":
        rows = [r for r in manifest if r["run_accession"] == args.run_accession]
        if len(rows) != 1:
            parser.error("trim requires one accession present in the manifest")
        trim_run(args.config, policy, scope, processing, rows[0])
    elif args.step == "finalize":
        finalize(args.config, policy, scope, manifest, processing)
    else:
        output = ROOT / policy["output_dir"]
        result = validate_artifacts(output, ROOT, [output / name for name in [*PRODUCTS, "config_snapshot.yaml", "input_checksums.json"]], policy)
        logging.info("Validated terminal QC: %s", result)


if __name__ == "__main__":
    main()
