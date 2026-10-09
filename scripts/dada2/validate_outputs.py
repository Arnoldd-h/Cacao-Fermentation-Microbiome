#!/usr/bin/env python3
"""Validate DADA2 pilot outputs using independent count and checksum checks.

Input: a completed run directory containing SUCCESS, checksum and ASV tables.
Output: JSON validation summary (stdout; optionally --output PATH).
Use --check-inputs-root to also verify the original inputs in a repository.
Fails on altered artifacts, missing/duplicate samples or ASVs, non-integer or
inconsistent counts, or a missing success marker. Does not interpret biology.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def read_rows(directory: Path, name: str) -> list[dict[str, str]]:
    with (directory / name).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty table: {name}")
    return rows


def integer(value: str) -> int:
    if not re.fullmatch(r"[0-9]+", value):
        raise ValueError(f"Expected a nonnegative integer count, got {value!r}")
    return int(value)


def unique_rows(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    result = {row[key]: row for row in rows}
    if len(result) != len(rows) or "" in result:
        raise ValueError(f"Missing or duplicate {key}")
    return result


def validate_run(directory: Path, inputs_root: Path | None = None) -> dict[str, object]:
    directory = directory.resolve()
    revision = (directory / "SUCCESS").read_text().strip()
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("Missing or invalid successful execution revision")
    checksums = read_rows(directory, "output_checksums.tsv")
    unique_rows(checksums, "path")
    for row in checksums:
        path = (directory / row["path"]).resolve()
        if directory not in path.parents:
            raise ValueError("Checksum output path escapes the run directory")
        with path.open("rb") as handle:
            observed = hashlib.file_digest(handle, "md5").hexdigest()
        if observed != row["md5"]:
            raise ValueError(f"Artifact checksum mismatch: {row['path']}")
    required = {"asv_counts.tsv", "asv_sequences.tsv", "asv_sequences.fasta", "read_tracking.tsv",
                "summary.tsv", "sample_metadata.tsv", "sequence_length_distribution.tsv", "error_learning.tsv",
                "input_checksums.tsv", "config_snapshot.yaml", "software_versions.tsv"}
    if not required.issubset({row["path"] for row in checksums}):
        raise ValueError("Checksum manifest is missing essential outputs")
    input_count = 0
    if inputs_root is not None:
        inputs_root = inputs_root.resolve()
        inputs = read_rows(directory, "input_checksums.tsv")
        unique_rows(inputs, "path")
        for row in inputs:
            path = (inputs_root / row["path"]).resolve()
            if inputs_root not in path.parents:
                raise ValueError("Checksum input path escapes the repository")
            with path.open("rb") as handle:
                observed = hashlib.file_digest(handle, "md5").hexdigest()
            if observed != row["md5"]:
                raise ValueError(f"Input checksum mismatch: {row['path']}")
        input_count = len(inputs)
    tracking = unique_rows(read_rows(directory, "read_tracking.tsv"), "sample_id")
    metadata = unique_rows(read_rows(directory, "sample_metadata.tsv"), "sample_id")
    counts = unique_rows(read_rows(directory, "asv_counts.tsv"), "sample_id")
    sequences = unique_rows(read_rows(directory, "asv_sequences.tsv"), "asv_id")
    if set(tracking) != set(metadata) or set(tracking) != set(counts):
        raise ValueError("Sample membership differs between metadata, tracking and ASV counts")
    studies = {row["study_id"] for table in (tracking, metadata, counts, sequences) for row in table.values()}
    if len(studies) != 1:
        raise ValueError("Outputs mix multiple studies")
    totals = {asv: 0 for asv in sequences}
    for sample_id, row in counts.items():
        if tracking[sample_id]["run_accession"] != metadata[sample_id]["run_accession"]:
            raise ValueError(f"Tracking run disagrees with metadata: {sample_id}")
        if set(row) - {"study_id", "sample_id"} != set(sequences):
            raise ValueError("ASV table columns disagree with sequence definitions")
        sample_total = 0
        for asv in sequences:
            abundance = integer(row[asv])
            totals[asv] += abundance
            sample_total += abundance
        if sample_total != integer(tracking[sample_id]["nonchim"]):
            raise ValueError(f"ASV counts disagree with retained pairs: {sample_id}")
        if sum(integer(row[asv]) > 0 for asv in sequences) != integer(tracking[sample_id]["asv_count"]):
            raise ValueError(f"ASV presence count disagrees with tracking: {sample_id}")
        steps = tracking[sample_id]
        values = {name: integer(steps[name]) for name in ("input", "filtered", "denoised_f", "denoised_r", "merged", "nonchim")}
        if not (0 < values["nonchim"] <= values["merged"] <= min(values["denoised_f"], values["denoised_r"])
                <= max(values["denoised_f"], values["denoised_r"]) <= values["filtered"] <= values["input"]):
            raise ValueError(f"Read conservation failed: {sample_id}")
    fasta = []
    for asv, row in sequences.items():
        if not re.fullmatch(r"[ACGT]+", row["sequence"]) or len(row["sequence"]) != integer(row["length"]):
            raise ValueError(f"Invalid ASV sequence or length: {asv}")
        if totals[asv] != integer(row["total_reads"]) or totals[asv] == 0:
            raise ValueError(f"ASV sequence total disagrees with counts: {asv}")
        fasta.extend([">" + asv, row["sequence"]])
    if (directory / "asv_sequences.fasta").read_text().splitlines() != fasta:
        raise ValueError("FASTA sequences disagree with the ASV table")
    summary = read_rows(directory, "summary.tsv")
    if len(summary) != 1 or integer(summary[0]["samples"]) != len(tracking) or integer(summary[0]["nonchim_asvs"]) != len(sequences):
        raise ValueError("Summary sample or ASV count mismatch")
    for column, stage in (("input_pairs", "input"), ("filtered_pairs", "filtered"), ("merged_pairs", "merged"), ("nonchim_pairs", "nonchim")):
        if integer(summary[0][column]) != sum(integer(row[stage]) for row in tracking.values()):
            raise ValueError(f"Summary count mismatch: {column}")
    lengths = read_rows(directory, "sequence_length_distribution.tsv")
    expected_lengths: dict[str, tuple[int, int]] = {}
    for row in sequences.values():
        asv_count, read_count = expected_lengths.get(row["length"], (0, 0))
        expected_lengths[row["length"]] = (asv_count + 1, read_count + integer(row["total_reads"]))
    retained_lengths = unique_rows([row for row in lengths if row["stage"] == "nonchim"], "length")
    observed_lengths = {length: (integer(row["asv_count"]), integer(row["read_count"])) for length, row in retained_lengths.items()}
    if observed_lengths != expected_lengths:
        raise ValueError("Non-chimeric length distribution disagrees with ASV sequences")
    merged_lengths = unique_rows([row for row in lengths if row["stage"] == "merged"], "length")
    if (sum(integer(row["asv_count"]) for row in merged_lengths.values()) != integer(summary[0]["merged_asvs"])
            or sum(integer(row["read_count"]) for row in merged_lengths.values()) != integer(summary[0]["merged_pairs"])):
        raise ValueError("Merged length distribution disagrees with the summary")
    errors = unique_rows(read_rows(directory, "error_learning.tsv"), "read_direction")
    if set(errors) != {"F", "R"}:
        raise ValueError("Expected forward and reverse error-learning diagnostics")
    return {"status": "valid", "source_git_commit": revision, "study_id": next(iter(studies)),
            "samples": len(tracking), "asvs": len(sequences), "nonchim_pairs": sum(totals.values()),
            "error_models_converged_before_limit": all(row["converged_before_limit"] == "TRUE" for row in errors.values()),
            "artifacts_checked": len(checksums), "inputs_checked": input_count}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=ROOT / "results/dada2/pilot")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-inputs-root", type=Path,
                        help="Also verify original input checksums relative to this repository root")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    summary = validate_run(args.run_dir, args.check_inputs_root)
    document = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(document, encoding="utf-8")
        temporary.replace(args.output)
    print(document, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
