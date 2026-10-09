"""Scoped QC paths and independently reconciled full-study read summaries."""

import json
from pathlib import Path

from .analysis_artifacts import read_table, unique_rows
from .pilot_workflow import fastq_path, processing_configuration
from .qc import build_cutadapt_summary_rows, reverse_complement_iupac
from .schema import CUTADAPT_SUMMARY_COLUMNS, PRIMER_DETECTION_COLUMNS


def trimmed_path(config: dict, row: dict, direction: str) -> str:
    return str(Path(config["interim_dir"]) / row["run_accession"] / (row["run_accession"] + "_" + direction + ".fastq.gz"))


def cutadapt_rows(config: dict, project: dict, manifest: list[dict], root: Path) -> list[dict]:
    rows = []
    for row in manifest:
        path = root / config["qc_dir"] / "cutadapt" / (row["run_accession"] + ".cutadapt.json")
        report = json.loads(path.read_text())
        processing = processing_configuration(project, row)
        outputs = {direction: Path(trimmed_path(config, row, suffix)).as_posix() for direction, suffix in (("R1", "1"), ("R2", "2"))}
        rows.extend(build_cutadapt_summary_rows(report, row, outputs, processing["primers"], processing["cutadapt"], root=root, report_path=path))
    return rows


def reconcile_qc(config: dict, project: dict, manifest: list[dict], root: Path) -> dict:
    directory = root / config["qc_dir"]
    expected = {(row["study_id"], row["run_accession"], direction): row for row in manifest for direction in ("R1", "R2")}
    fields = ("study_id", "run_accession", "read_direction")
    raw = unique_rows(read_table(directory / "raw_read_quality.tsv", list(fields) + ["total_sequences", "sample_id"]), fields)
    trimmed = unique_rows(read_table(directory / "trimmed_read_quality.tsv", list(fields) + ["total_sequences", "sample_id"]), fields)
    cutadapt = unique_rows(read_table(directory / "cutadapt_summary.tsv", CUTADAPT_SUMMARY_COLUMNS), fields)
    independently_parsed = unique_rows(cutadapt_rows(config, project, manifest, root), fields)
    if set(raw) != set(expected) or set(trimmed) != set(expected) or cutadapt != independently_parsed:
        raise ValueError("QC run membership or Cutadapt parameters/counts changed")
    total = kept = 0
    for key, sample in expected.items():
        for table in (raw, trimmed, cutadapt):
            if table[key]["sample_id"] != sample["sample_id"]:
                raise ValueError("QC sample metadata disagrees with full manifest")
        if int(raw[key]["total_sequences"]) != int(cutadapt[key]["input_reads"]) or int(trimmed[key]["total_sequences"]) != int(cutadapt[key]["reads_written"]):
            raise ValueError("Raw/trimmed FastQC and Cutadapt counts disagree")
        other = (key[0], key[1], "R2" if key[2] == "R1" else "R1")
        if raw[key]["total_sequences"] != raw[other]["total_sequences"] or trimmed[key]["total_sequences"] != trimmed[other]["total_sequences"]:
            raise ValueError("Unequal paired read counts")
        if key[2] == "R1":
            total += int(raw[key]["total_sequences"])
            kept += int(trimmed[key]["total_sequences"])
    for stage, quality in (("raw", raw), ("trimmed", trimmed)):
        screening = read_table(directory / ("primer_detection_" + stage + ".tsv"), PRIMER_DETECTION_COLUMNS)
        keys = ("study_id", "run_accession", "read_direction", "primer", "sequence_type", "orientation")
        observed = unique_rows(screening, keys)
        required = {}
        for key, sample in expected.items():
            processing = processing_configuration(project, sample)
            for role, primer in processing["primers"].items():
                label = primer["name"].split("_", 1)[0]
                for kind, field in (("gene_specific", "gene_specific_sequence"), ("full_construct", "full_construct_sequence")):
                    sequence = primer[field]
                    for orientation in ("as_synthesized", "reverse_complement"):
                        required[(*key, label, kind, orientation)] = {
                            "primer_sequence": sequence,
                            "search_sequence": sequence if orientation == "as_synthesized" else reverse_complement_iupac(sequence),
                            "error_rate": str(processing["primer_detection"]["error_rate"]),
                            "minimum_overlap": str(len(sequence)),
                            "search_scope": "five_prime_anchored" if kind == "full_construct" else "any_position_full_length",
                            "expected_case": str(orientation == "as_synthesized" and ((role == "forward" and key[2] == "R1") or (role == "reverse" and key[2] == "R2"))).lower()}
        if set(observed) != set(required):
            raise ValueError("Incomplete or unexpected primer screening")
        for row in screening:
            key = tuple(row[field] for field in fields)
            if key not in expected or int(row["reads_examined"]) != int(quality[key]["total_sequences"]) or not 0 <= int(row["primer_matches"]) <= int(row["reads_examined"]):
                raise ValueError("Primer screen counts disagree with QC")
            if row["sample_id"] != expected[key]["sample_id"] or any(row[field] != value for field, value in required[tuple(row[field] for field in keys)].items()):
                raise ValueError("Primer screen metadata or executed parameters disagree")
            if abs(float(row["match_percent"]) - int(row["primer_matches"]) / int(row["reads_examined"]) * 100) > 1e-6:
                raise ValueError("Primer match percentage disagrees with counts")
    return {"study_id": config["study_id"], "runs": len(manifest), "fastq_files": len(expected), "input_pairs": total,
            "trimmed_pairs": kept, "paired_retention_fraction": kept / total, "sample_exclusions": 0,
            "raw_base_quality_pass_files": sum(row["per_base_quality_status"] == "pass" for row in raw.values()),
            "trimmed_base_quality_pass_files": sum(row["per_base_quality_status"] == "pass" for row in trimmed.values())}
