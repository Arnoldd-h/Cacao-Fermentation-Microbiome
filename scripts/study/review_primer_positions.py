#!/usr/bin/env python3
"""Inspect the largest opposite-primer RC residual case in one trimmed mate.

Inputs: scope/profile, validated manifest, trimmed primer screen and FASTQ.
Outputs: Cutadapt JSON, aggregate match positions, summary and provenance.
Intermediate per-read info stays under results/intermediate (not versioned).
Fails on invalid coordinates or counts inconsistent with the original screen.
This is a selected illustrative QC case, not a cohort-wide estimate or filter.
"""

import argparse
from collections import Counter
import json
import logging
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.analysis_artifacts import read_table, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.pilot_workflow import load_pilot_manifest, manifest_row, processing_configuration
from cacao_inventory.provenance import write_provenance
from cacao_inventory.qc import cutadapt_detection_counts, reverse_complement_iupac
from cacao_inventory.study_qc import trimmed_path
from cacao_inventory.study_scope import validate_scope


def read_positions(path):
    total = 0
    positions = Counter()
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            total += 1
            if len(fields) >= 4 and fields[1] == "-1":
                if len(fields[2]) != len(fields[3]):
                    raise ValueError("Unmatched info sequence/quality lengths disagree")
                continue
            if len(fields) < 8 or int(fields[1]) < 0:
                raise ValueError("Malformed matched Cutadapt info row")
            start, end = int(fields[2]), int(fields[3])
            length = len(fields[4] + fields[5] + fields[6])
            if start != len(fields[4]) or end - start != len(fields[5]) or not 0 <= start < end <= length:
                raise ValueError("Cutadapt info coordinates disagree with sequence segments")
            positions[(start, end, length)] += 1
    return total, positions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--direction", choices=("R1", "R2"), required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    scope = load_json_yaml(args.config)
    validate_scope(scope, ROOT)
    qc_validation = ROOT / scope["qc_dir"] / "validation.json"
    checked = load_json_yaml(qc_validation)
    if checked["status"] != "valid" or checked["study_id"] != scope["study_id"]:
        raise ValueError("A completed QC validation for this study is required")
    project = load_json_yaml(ROOT / scope["project_config"])
    manifest = load_pilot_manifest(ROOT / scope["manifest"])
    processing = processing_configuration(project, manifest[0])
    role = "reverse" if args.direction == "R1" else "forward"
    label = processing["primers"][role]["name"].split("_", 1)[0]
    screen = ROOT / scope["qc_dir"] / "primer_detection_trimmed.tsv"
    cases = [row for row in read_table(screen, []) if row["read_direction"] == args.direction
             and row["primer"] == label and row["orientation"] == "reverse_complement"
             and row["sequence_type"] == "gene_specific"]
    if not cases:
        raise ValueError("No opposite-primer RC cases in the completed screen")
    case = sorted(cases, key=lambda row: (-int(row["primer_matches"]), row["run_accession"]))[0]
    sample = manifest_row(manifest, case["run_accession"], scope["study_id"])
    read = ROOT / trimmed_path(scope, sample, "1" if args.direction == "R1" else "2")
    output = ROOT / scope["qc_dir"] / "position_review" / args.direction
    intermediate = ROOT / "results/intermediate/qc" / scope["study_id"] / "position_review"
    output.mkdir(parents=True, exist_ok=True)
    intermediate.mkdir(parents=True, exist_ok=True)
    info = intermediate / (args.direction + ".info.tsv")
    report = output / "cutadapt.json"
    sequence = reverse_complement_iupac(processing["primers"][role]["gene_specific_sequence"])
    command = ["cutadapt", "--action=none", "--no-indels", "--error-rate", str(processing["primer_detection"]["error_rate"]),
               "--overlap", str(len(sequence)), "--front", "probe=" + sequence,
               "--info-file", str(info.relative_to(ROOT)), "--json", str(report.relative_to(ROOT)),
               "--output", os.devnull, str(read.relative_to(ROOT))]
    subprocess.run(command, cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    result = json.loads(report.read_text())
    examined, matched = cutadapt_detection_counts(result)
    total, positions = read_positions(info)
    if (examined, matched) != (int(case["reads_examined"]), int(case["primer_matches"])) or total != examined or sum(positions.values()) != matched:
        raise ValueError("Position replay disagrees with the original primer screen")
    rows = [{"match_start_0_based": start, "match_end_exclusive": end, "read_length": length,
             "bases_after_match": length - end, "read_count": count} for (start, end, length), count in sorted(positions.items())]
    columns = ["match_start_0_based", "match_end_exclusive", "read_length", "bases_after_match", "read_count"]
    write_tsv_atomic(output / "positions.tsv", rows, columns)
    terminal = sum(count for (start, end, length), count in positions.items() if end == length)
    summary = {"study_id": scope["study_id"], "case": case,
               "selection": "largest matched count for the opposite primer RC in this mate; ties by run accession",
               "illustrative_not_cohort_estimate": True, "coordinate_system": "zero_based_half_open",
               "reads_examined": examined, "matched_reads": matched, "terminal_matches": terminal,
               "terminal_fraction_of_matches": terminal / matched if matched else None,
               "read_filtering": False, "sample_exclusions": 0}
    write_json(output / "summary.json", summary)
    write_provenance(output / "provenance.json", root=ROOT,
                     inputs=[screen, qc_validation, ROOT / scope["qc_dir"] / "provenance.json",
                             ROOT / scope["manifest"], read, Path(__file__), ROOT / "environment/conda-linux-64.lock",
                             *[ROOT / "python/cacao_inventory" / (name + ".py") for name in
                               ("qc", "config", "io", "provenance", "pilot_workflow", "study_scope", "study_qc", "analysis_artifacts")]],
                     outputs=[report, output / "positions.tsv", output / "summary.json"],
                     config_paths=[args.config, ROOT / scope["project_config"]],
                     software_versions={"Cutadapt": result["cutadapt_version"]}, command=sys.argv)
    logging.info("%s %s: %d/%d matches are terminal", args.direction, case["run_accession"], terminal, matched)


if __name__ == "__main__":
    main()
