"""Quality-control transformations for amplicon sequencing reports."""

from __future__ import annotations

from decimal import Decimal, InvalidOperation
from pathlib import Path


MULTIQC_FASTQC_COLUMNS = [
    "Sample",
    "Filename",
    "Total Sequences",
    "Sequences flagged as poor quality",
    "Sequence length",
    "%GC",
    "avg_sequence_length",
    "median_sequence_length",
    "basic_statistics",
    "per_base_sequence_quality",
    "per_sequence_quality_scores",
    "per_base_sequence_content",
    "per_sequence_gc_content",
    "per_base_n_content",
    "sequence_length_distribution",
    "sequence_duplication_levels",
    "overrepresented_sequences",
    "adapter_content",
]


def _integer(value: str, field: str) -> int:
    """Parse integer-like values, including the decimal form emitted by MultiQC."""

    try:
        number = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError(f"{field} is not numeric: {value}") from exc
    if number != number.to_integral_value():
        raise ValueError(f"{field} is not an integer: {value}")
    return int(number)


def read_direction(read_file: str) -> str:
    """Return R1 or R2 from a paired ENA FASTQ filename."""

    name = Path(read_file).name
    if name.endswith("_1.fastq.gz"):
        return "R1"
    if name.endswith("_2.fastq.gz"):
        return "R2"
    raise ValueError(f"Cannot determine paired read direction: {read_file}")


def fastqc_sample_id(read_file: str) -> str:
    """Return the sample identifier used by FastQC and MultiQC."""

    suffix = ".fastq.gz"
    name = Path(read_file).name
    if not name.endswith(suffix):
        raise ValueError(f"Expected a .fastq.gz filename: {read_file}")
    return name[: -len(suffix)]


def build_raw_quality_rows(
    multiqc_rows: list[dict[str, str]],
    validation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Join MultiQC FastQC metrics to validated FASTQ provenance."""

    fastqc_by_sample: dict[str, dict[str, str]] = {}
    for row in multiqc_rows:
        sample = row["Sample"]
        if sample in fastqc_by_sample:
            raise ValueError(f"Duplicate MultiQC FastQC sample: {sample}")
        fastqc_by_sample[sample] = row

    validation_by_sample: dict[str, dict[str, str]] = {}
    for row in validation_rows:
        sample = fastqc_sample_id(row["read_file"])
        if sample in validation_by_sample:
            raise ValueError(f"Duplicate validated FASTQ sample: {sample}")
        validation_by_sample[sample] = row

    missing = sorted(set(validation_by_sample) - set(fastqc_by_sample))
    unexpected = sorted(set(fastqc_by_sample) - set(validation_by_sample))
    if missing or unexpected:
        raise ValueError(
            "FastQC/FASTQ sample mismatch; "
            f"missing={','.join(missing) or 'none'}; "
            f"unexpected={','.join(unexpected) or 'none'}"
        )

    output: list[dict[str, str]] = []
    for sample in sorted(validation_by_sample):
        validation = validation_by_sample[sample]
        fastqc = fastqc_by_sample[sample]
        if validation["status"] != "valid":
            raise ValueError(f"FASTQ was not validated before FastQC: {sample}")
        if fastqc["Filename"] != validation["read_file"]:
            raise ValueError(f"FastQC filename mismatch for {sample}")
        total_sequences = _integer(fastqc["Total Sequences"], "Total Sequences")
        expected_sequences = _integer(validation["read_count"], "read_count")
        if total_sequences != expected_sequences:
            raise ValueError(
                f"FastQC read count mismatch for {sample}: "
                f"{total_sequences} != {expected_sequences}"
            )
        output.append(
            {
                "study_id": validation["study_id"],
                "bioproject": validation["bioproject"],
                "sample_id": validation["sample_id"],
                "run_accession": validation["run_accession"],
                "read_direction": read_direction(validation["read_file"]),
                "read_file": validation["read_file"],
                "total_sequences": str(total_sequences),
                "total_bases": validation["base_count"],
                "minimum_read_length": validation["minimum_read_length"],
                "maximum_read_length": validation["maximum_read_length"],
                "sequence_length": fastqc["Sequence length"],
                "average_sequence_length": fastqc["avg_sequence_length"],
                "median_sequence_length": fastqc["median_sequence_length"],
                "gc_percent": fastqc["%GC"],
                "poor_quality_sequences": str(
                    _integer(
                        fastqc["Sequences flagged as poor quality"],
                        "Sequences flagged as poor quality",
                    )
                ),
                "basic_statistics_status": fastqc["basic_statistics"],
                "per_base_quality_status": fastqc["per_base_sequence_quality"],
                "per_sequence_quality_status": fastqc["per_sequence_quality_scores"],
                "per_base_content_status": fastqc["per_base_sequence_content"],
                "per_sequence_gc_status": fastqc["per_sequence_gc_content"],
                "per_base_n_content_status": fastqc["per_base_n_content"],
                "sequence_length_distribution_status": fastqc[
                    "sequence_length_distribution"
                ],
                "sequence_duplication_status": fastqc["sequence_duplication_levels"],
                "overrepresented_sequences_status": fastqc[
                    "overrepresented_sequences"
                ],
                "adapter_content_status": fastqc["adapter_content"],
            }
        )

    rows_by_run: dict[str, list[dict[str, str]]] = {}
    for row in output:
        rows_by_run.setdefault(row["run_accession"], []).append(row)
    for run_accession, run_rows in rows_by_run.items():
        directions = {row["read_direction"] for row in run_rows}
        if directions != {"R1", "R2"}:
            raise ValueError(
                f"Paired FASTQ directions are incomplete for {run_accession}: "
                f"{','.join(sorted(directions)) or 'none'}"
            )
        counts = {row["total_sequences"] for row in run_rows}
        if len(counts) != 1:
            raise ValueError(f"Paired read counts disagree for {run_accession}")
        sample_ids = {row["sample_id"] for row in run_rows}
        if len(sample_ids) != 1:
            raise ValueError(f"Paired sample IDs disagree for {run_accession}")
    return output
