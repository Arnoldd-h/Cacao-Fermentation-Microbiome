#!/usr/bin/env python3
"""Execute one scoped QC step using the registered, unchanged primer methods.

Inputs: full-study config/manifest and declared FASTQ/Cutadapt/FastQC inputs.
Outputs: per-run trimming/screens, aggregate tables or QC hash/validation record.
Fails on altered scope, unsupported methods, tool errors or inconsistent counts.
"""

import argparse
import json
import logging
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/qc"))
from detect_primers import _processing_config, _run_cutadapt
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, read_table, validate_artifacts, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.pilot_workflow import fastq_path, load_pilot_manifest, manifest_row, processing_configuration
from cacao_inventory.provenance import collect_versions, file_records
from cacao_inventory.qc import cutadapt_command, reverse_complement_iupac, validate_iupac_sequence
from cacao_inventory.schema import CUTADAPT_SUMMARY_COLUMNS, PRIMER_DETECTION_COLUMNS
from cacao_inventory.study_qc import cutadapt_rows, reconcile_qc, trimmed_path
from cacao_inventory.study_scope import validate_scope

PRODUCTS = ["download_validation.tsv", "fastq_validation.tsv", "raw_read_quality.tsv", "trimmed_read_quality.tsv",
            "read_quality_comparison.tsv", "cutadapt_summary.tsv", "primer_detection_raw.tsv", "primer_detection_trimmed.tsv", "summary.json"]


def screen(config, project, row, stage):
    processing = _processing_config(project, row["bioproject"])
    rows = []
    for direction, suffix in (("R1", "1"), ("R2", "2")):
        path = ROOT / (fastq_path(row, suffix) if stage == "raw" else trimmed_path(config, row, suffix))
        for role, primer in processing["primers"].items():
            label = primer["name"].split("_", 1)[0]
            for kind, field, anchored in (("gene_specific", "gene_specific_sequence", False), ("full_construct", "full_construct_sequence", True)):
                sequence = validate_iupac_sequence(primer[field], field)
                for orientation, search in (("as_synthesized", sequence), ("reverse_complement", reverse_complement_iupac(sequence))):
                    examined, matched = _run_cutadapt("cutadapt", path, label + "_" + kind + "_" + orientation, search, anchored, processing["primer_detection"]["error_rate"])
                    rows.append({"study_id": row["study_id"], "bioproject": row["bioproject"], "sample_id": row["sample_id"], "run_accession": row["run_accession"],
                        "read_direction": direction, "read_file": path.name, "primer": label, "sequence_type": kind, "primer_sequence": sequence,
                        "search_sequence": search, "reads_examined": str(examined), "primer_matches": str(matched), "match_percent": f"{matched / examined * 100:.6f}",
                        "orientation": orientation, "search_scope": "five_prime_anchored" if anchored else "any_position_full_length",
                        "error_rate": str(processing["primer_detection"]["error_rate"]), "minimum_overlap": str(len(sequence)),
                        "expected_case": str(orientation == "as_synthesized" and ((role == "forward" and direction == "R1") or (role == "reverse" and direction == "R2"))).lower()})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("trim", "screen", "summarize_cutadapt", "summarize_screen", "finalize"))
    parser.add_argument("--config", type=Path, default=ROOT / "config/full_study.yaml")
    parser.add_argument("--run-accession")
    parser.add_argument("--stage", choices=("raw", "trimmed"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    config = load_json_yaml(args.config)
    validate_scope(config, ROOT)
    project = load_json_yaml(ROOT / config["project_config"])
    manifest = load_pilot_manifest(ROOT / config["manifest"])
    directory = ROOT / config["qc_dir"]
    if args.step in ("trim", "screen"):
        if not args.run_accession:
            parser.error("Per-run steps require --run-accession")
        row = manifest_row(manifest, args.run_accession, config["study_id"])
        if args.step == "trim":
            processing = processing_configuration(project, row)
            outputs = {direction: Path(trimmed_path(config, row, suffix)).as_posix() for direction, suffix in (("R1", "1"), ("R2", "2"))}
            report = directory / "cutadapt" / (args.run_accession + ".cutadapt.json")
            for path in [*[ROOT / value for value in outputs.values()], report]:
                path.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(cutadapt_command(row, processing["primers"], processing["cutadapt"], outputs, str(report.relative_to(ROOT)), 1), cwd=ROOT, check=True)
        else:
            if not args.stage or not args.output:
                parser.error("Screening requires --stage and --output")
            write_tsv_atomic(args.output, screen(config, project, row, args.stage), PRIMER_DETECTION_COLUMNS)
    elif args.step == "summarize_cutadapt":
        write_tsv_atomic(directory / "cutadapt_summary.tsv", cutadapt_rows(config, project, manifest, ROOT), CUTADAPT_SUMMARY_COLUMNS)
    elif args.step == "summarize_screen":
        if not args.stage:
            parser.error("Screen aggregation requires --stage")
        rows = []
        for row in manifest:
            path = ROOT / "results/intermediate/qc" / config["study_id"] / ("primers_" + args.stage) / (row["run_accession"] + ".tsv")
            rows.extend(read_table(path, PRIMER_DETECTION_COLUMNS))
        write_tsv_atomic(directory / ("primer_detection_" + args.stage + ".tsv"), rows, PRIMER_DETECTION_COLUMNS)
    else:
        state = begin_stage(directory, ROOT)
        inputs = [args.config, ROOT / config["project_config"], ROOT / config["manifest"], Path(__file__),
                  ROOT / "python/cacao_inventory/study_qc.py", ROOT / "workflow/study.Snakefile", ROOT / "environment/conda-linux-64.lock", *[directory / name for name in PRODUCTS if name != "summary.json"]]
        inputs += list((ROOT / "python/cacao_inventory").glob("*.py")) + [ROOT / "scripts/qc/detect_primers.py", ROOT / "scripts/qc/summarize_fastqc.py", ROOT / "scripts/qc/summarize_trimmed_fastqc.py"]
        for row in manifest:
            inputs += [directory / "cutadapt" / (row["run_accession"] + ".cutadapt.json")]
            inputs += [ROOT / fastq_path(row, d) for d in ("1", "2")]
            inputs += [ROOT / trimmed_path(config, row, d) for d in ("1", "2")]
        before = file_records(inputs, ROOT)
        summary = reconcile_qc(config, project, manifest, ROOT)
        write_json(directory / "summary.json", summary)
        products = [directory / name for name in PRODUCTS]
        finish_stage(directory, ROOT, state, before, args.config, config, products, {}, collect_versions("fastqc", "multiqc", "cutadapt"), sys.argv)
        validation = validate_artifacts(directory, ROOT, [*products, directory / "config_snapshot.yaml", directory / "input_checksums.json"], config)
        write_json(directory / "validation.json", {**validation, **summary})
        logging.info("Validated full-study QC: %s", summary)


if __name__ == "__main__":
    main()
