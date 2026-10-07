"""Quality-control transformations for amplicon sequencing reports."""

from __future__ import annotations

from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any
import math

from .pilot_workflow import fastq_path


IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)


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


def validate_iupac_sequence(sequence: str, field: str = "sequence") -> str:
    """Return an uppercase IUPAC DNA sequence or fail with a useful message."""

    normalized = sequence.upper()
    if not normalized:
        raise ValueError(f"{field} cannot be empty")
    invalid = sorted(set(normalized) - IUPAC_DNA)
    if invalid:
        raise ValueError(f"{field} contains invalid IUPAC symbols: {''.join(invalid)}")
    return normalized


def reverse_complement_iupac(sequence: str) -> str:
    """Reverse-complement an IUPAC DNA sequence, preserving ambiguity codes."""

    normalized = validate_iupac_sequence(sequence)
    return normalized.translate(IUPAC_COMPLEMENT)[::-1]


def cutadapt_detection_counts(report: dict[str, Any]) -> tuple[int, int]:
    """Extract and cross-check the examined and adapter-matched read counts."""

    try:
        read_counts = report["read_counts"]
        examined = int(read_counts["input"])
        matched = int(read_counts["read1_with_adapter"])
        adapters = report["adapters_read1"]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Cutadapt JSON lacks primer-detection counts") from exc
    if not isinstance(adapters, list) or len(adapters) != 1:
        raise ValueError("Primer detection requires exactly one Cutadapt adapter")
    try:
        adapter_matches = int(adapters[0]["total_matches"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Cutadapt JSON lacks adapter total_matches") from exc
    if examined < 0 or matched < 0 or matched > examined:
        raise ValueError("Cutadapt primer-detection read counts are inconsistent")
    if adapter_matches != matched:
        raise ValueError(
            "Cutadapt adapter/read match counts disagree: "
            f"{adapter_matches} != {matched}"
        )
    return examined, matched


def validate_cutadapt_parameters(primers: dict[str, Any], parameters: dict[str, Any]) -> None:
    """Reject unsupported method changes instead of silently ignoring flags."""

    for flag in ("allow_indels", "discard_untrimmed", "quality_trimming"):
        if parameters.get(flag) is not False:
            raise ValueError(f"Pilot Cutadapt requires {flag}=false; other values are unsupported")
    error_rate = float(parameters["error_rate"])
    if not math.isfinite(error_rate) or not 0 <= error_rate <= 1:
        raise ValueError("Invalid Cutadapt error_rate")
    if set(primers) != {"forward", "reverse"}:
        raise ValueError("Cutadapt requires forward and reverse primers")
    for role in ("forward", "reverse"):
        sequence = validate_iupac_sequence(primers[role]["gene_specific_sequence"])
        if int(parameters[f"{role}_minimum_overlap"]) != len(sequence):
            raise ValueError(f"Cutadapt requires full-length {role} primer overlap")


def cutadapt_command(
    row: dict[str, str], primers: dict[str, Any], parameters: dict[str, Any],
    outputs: dict[str, str], report_path: str, threads: int = 1,
) -> list[str]:
    """Construct the only currently supported primer-only paired trimming command."""

    validate_cutadapt_parameters(primers, parameters)
    if threads < 1:
        raise ValueError("Cutadapt threads must be positive")
    adapters = {
        role: f"{role}={validate_iupac_sequence(primers[role]['gene_specific_sequence'])};min_overlap={parameters[f'{role}_minimum_overlap']}"
        for role in ("forward", "reverse")
    }
    return [
        "cutadapt", "--cores", str(threads), "--no-indels", "--error-rate", str(parameters["error_rate"]),
        "--front", adapters["forward"], "-G", adapters["reverse"],
        "--json", report_path, "--output", outputs["R1"], "--paired-output", outputs["R2"],
        fastq_path(row, "1"), fastq_path(row, "2"),
    ]


def _validate_cutadapt_execution(
    report: dict[str, Any], row: dict[str, str], outputs: dict[str, str],
    primers: dict[str, Any], parameters: dict[str, Any], root: Path,
    report_path: str | Path | None,
) -> None:
    """Compare actual report paths, command flags and adapters to configured inputs."""

    validate_cutadapt_parameters(primers, parameters)
    arguments = report.get("command_line_arguments")
    if not isinstance(arguments, list):
        raise ValueError("Cutadapt JSON lacks command_line_arguments")
    aliases = {"-j": "--cores", "-g": "--front", "-o": "--output", "-p": "--paired-output"}
    valued = {"--cores", "--error-rate", "--front", "-G", "--json", "--output", "--paired-output"}
    options: dict[str, str | bool] = {}
    positional: list[str] = []
    cursor = 0
    while cursor < len(arguments):
        token = aliases.get(arguments[cursor], arguments[cursor])
        if token in options:
            raise ValueError(f"Duplicate Cutadapt option: {token}")
        if token == "--no-indels":
            options[token] = True
        elif token in valued:
            cursor += 1
            if cursor >= len(arguments):
                raise ValueError(f"Missing Cutadapt option value: {token}")
            options[token] = str(arguments[cursor])
        elif token.startswith("-"):
            raise ValueError(f"Unsupported Cutadapt command option: {token}")
        else:
            positional.append(token)
        cursor += 1
    if set(options) != valued | {"--no-indels"} or len(positional) != 2:
        raise ValueError("Cutadapt command does not match the paired primer-only method")
    if int(options["--cores"]) < 1 or float(options["--error-rate"]) != float(parameters["error_rate"]):
        raise ValueError("Cutadapt command error rate/cores disagrees with configuration")

    def canonical(value: str | Path) -> Path:
        return (root / value).resolve()

    for suffix, direction, path_field, output_option, role, adapter_option, adapter_key in (
        ("1", "R1", "path1", "--output", "forward", "--front", "adapters_read1"),
        ("2", "R2", "path2", "--paired-output", "reverse", "-G", "adapters_read2"),
    ):
        expected_input = canonical(fastq_path(row, suffix))
        expected_output = canonical(fastq_path(row, suffix, "trimmed"))
        if canonical(report["input"][path_field]) != expected_input or canonical(positional[int(suffix) - 1]) != expected_input:
            raise ValueError(f"Cutadapt input path disagrees with manifest for {direction}")
        if canonical(outputs[direction]) != expected_output or canonical(str(options[output_option])) != expected_output:
            raise ValueError(f"Cutadapt output path disagrees with manifest for {direction}")
        sequence = validate_iupac_sequence(primers[role]["gene_specific_sequence"])
        adapter_value = str(options[adapter_option])
        name, separator, specification = adapter_value.partition("=")
        expected_specification = f"{sequence};min_overlap={parameters[f'{role}_minimum_overlap']}"
        if not separator or not name or specification != expected_specification:
            raise ValueError(f"Cutadapt command primer/overlap disagrees for {direction}")
        adapters = report.get(adapter_key)
        if not isinstance(adapters, list) or len(adapters) != 1:
            raise ValueError(f"Cutadapt JSON must contain one adapter for {direction}")
        adapter = adapters[0]
        end = adapter.get("five_prime_end") or {}
        if (adapter.get("name") != name or end.get("sequence") != sequence
                or end.get("indels") is not False
                or end.get("type") != "regular_five_prime"
                or float(end.get("error_rate", -1)) != float(parameters["error_rate"])):
            raise ValueError(f"Cutadapt reported primer/parameters disagree for {direction}")
        count_key = "read1_with_adapter" if suffix == "1" else "read2_with_adapter"
        if int(adapter["total_matches"]) != int(report["read_counts"][count_key]):
            raise ValueError(f"Cutadapt reported adapter counts disagree for {direction}")
    if report_path is not None and canonical(str(options["--json"])) != canonical(report_path):
        raise ValueError("Cutadapt JSON output path disagrees with report location")
    if any(value is not None for value in report["read_counts"].get("filtered", {}).values()):
        raise ValueError("Cutadapt report unexpectedly enabled read filtering")
    if report["basepair_counts"].get("quality_trimmed") is not None:
        raise ValueError("Cutadapt report unexpectedly enabled quality trimming")


def build_cutadapt_summary_rows(
    report: dict[str, Any],
    manifest_row: dict[str, str],
    output_files: dict[str, str],
    primers: dict[str, dict[str, Any]],
    parameters: dict[str, Any],
    *,
    root: Path | None = None,
    report_path: str | Path | None = None,
) -> list[dict[str, str]]:
    """Validate a paired Cutadapt JSON report and return one row per direction."""

    try:
        input_report = report["input"]
        read_counts = report["read_counts"]
        base_counts = report["basepair_counts"]
        version = str(report["cutadapt_version"])
        input_pairs = int(read_counts["input"])
        output_pairs = int(read_counts["output"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Cutadapt JSON lacks paired trimming metrics") from exc
    if input_report.get("paired") is not True:
        raise ValueError("Cutadapt pilot trimming report is not paired-end")
    if input_pairs <= 0 or output_pairs < 0 or output_pairs > input_pairs:
        raise ValueError("Cutadapt paired read counts are inconsistent")
    if parameters.get("discard_untrimmed") is not False:
        raise ValueError("Pilot trimming must retain reads without verified primers")
    if parameters.get("allow_indels") is not False:
        raise ValueError("Pilot primer trimming must disable indels")
    if parameters.get("quality_trimming") is not False:
        raise ValueError("Quality trimming is not part of the Cutadapt primer step")
    if output_pairs != input_pairs:
        raise ValueError(
            "Cutadapt unexpectedly removed paired reads despite non-discarding configuration"
        )

    _validate_cutadapt_execution(
        report, manifest_row, output_files, primers, parameters, root or Path.cwd(), report_path,
    )

    try:
        error_rate = float(parameters["error_rate"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Invalid Cutadapt error rate in configuration") from exc

    rows: list[dict[str, str]] = []
    direction_fields = {
        "R1": (
            "path1",
            "input_read1",
            "read1_with_adapter",
            "output_read1",
            "forward",
            "forward_minimum_overlap",
        ),
        "R2": (
            "path2",
            "input_read2",
            "read2_with_adapter",
            "output_read2",
            "reverse",
            "reverse_minimum_overlap",
        ),
    }
    for direction, fields in direction_fields.items():
        path_field, input_bases_field, matches_field, output_bases_field, role, overlap_field = (
            fields
        )
        try:
            input_file = str(input_report[path_field])
            input_bases = int(base_counts[input_bases_field])
            matches = int(read_counts[matches_field])
            output_bases = int(base_counts[output_bases_field])
            primer = primers[role]
            primer_sequence = validate_iupac_sequence(
                str(primer["gene_specific_sequence"]), f"{role} primer"
            )
            minimum_overlap = int(parameters[overlap_field])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Cutadapt JSON/config lacks {direction} metrics") from exc
        if not 0 <= matches <= input_pairs:
            raise ValueError(f"Cutadapt adapter matches are inconsistent for {direction}")
        if input_bases <= 0 or not 0 <= output_bases <= input_bases:
            raise ValueError(f"Cutadapt base counts are inconsistent for {direction}")
        if minimum_overlap != len(primer_sequence):
            raise ValueError(f"Cutadapt must require the full {direction} primer sequence")
        rows.append(
            {
                "study_id": manifest_row["study_id"],
                "bioproject": manifest_row["bioproject"],
                "sample_id": manifest_row["sample_id"],
                "run_accession": manifest_row["run_accession"],
                "read_direction": direction,
                "input_file": input_file,
                "output_file": output_files[direction],
                "input_reads": str(input_pairs),
                "input_bases": str(input_bases),
                "reads_with_adapter": str(matches),
                "reads_written": str(output_pairs),
                "output_bases": str(output_bases),
                "percent_retained": f"{output_pairs / input_pairs * 100:.6f}",
                "primer": str(primer["name"]),
                "primer_sequence": primer_sequence,
                "error_rate": str(error_rate),
                "minimum_overlap": str(minimum_overlap),
                "discard_untrimmed": "false",
                "cutadapt_version": version,
            }
        )
    return rows


def fastqc_length_bounds(value: str) -> tuple[int, int]:
    """Parse FastQC's integer or inclusive range sequence-length field."""

    parts = value.split("-", maxsplit=1)
    try:
        lower = int(parts[0])
        upper = int(parts[-1])
    except ValueError as exc:
        raise ValueError(f"Invalid FastQC sequence length: {value}") from exc
    if lower <= 0 or upper < lower:
        raise ValueError(f"Invalid FastQC sequence length: {value}")
    return lower, upper


def build_trimmed_quality_rows(
    multiqc_rows: list[dict[str, str]],
    cutadapt_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Join trimmed FastQC metrics to validated Cutadapt counts and provenance."""

    fastqc_by_sample = {row["Sample"]: row for row in multiqc_rows}
    if len(fastqc_by_sample) != len(multiqc_rows):
        raise ValueError("Duplicate trimmed MultiQC FastQC sample")
    pseudo_validation: list[dict[str, str]] = []
    output_path_by_key: dict[tuple[str, str], str] = {}
    for row in cutadapt_rows:
        output_name = Path(row["output_file"]).name
        sample = fastqc_sample_id(output_name)
        if sample not in fastqc_by_sample:
            raise ValueError(f"Missing trimmed FastQC sample: {sample}")
        lower, upper = fastqc_length_bounds(fastqc_by_sample[sample]["Sequence length"])
        direction = row["read_direction"]
        key = (row["run_accession"], direction)
        if key in output_path_by_key:
            raise ValueError(f"Duplicate Cutadapt summary row: {key}")
        output_path_by_key[key] = row["output_file"]
        pseudo_validation.append(
            {
                "study_id": row["study_id"],
                "bioproject": row["bioproject"],
                "sample_id": row["sample_id"],
                "run_accession": row["run_accession"],
                "read_file": output_name,
                "read_count": row["reads_written"],
                "base_count": row["output_bases"],
                "minimum_read_length": str(lower),
                "maximum_read_length": str(upper),
                "status": "valid",
            }
        )
    output = build_raw_quality_rows(multiqc_rows, pseudo_validation)
    for row in output:
        key = (row["run_accession"], row["read_direction"])
        row["read_file"] = output_path_by_key[key]
    return output


def build_read_quality_comparison_rows(
    raw_rows: list[dict[str, str]],
    trimmed_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Compare raw and post-Cutadapt FastQC metrics by run and direction."""

    def index(rows: list[dict[str, str]], label: str) -> dict[tuple[str, str], dict[str, str]]:
        indexed: dict[tuple[str, str], dict[str, str]] = {}
        for row in rows:
            key = (row["run_accession"], row["read_direction"])
            if key in indexed:
                raise ValueError(f"Duplicate {label} QC row: {key}")
            indexed[key] = row
        return indexed

    raw_by_key = index(raw_rows, "raw")
    trimmed_by_key = index(trimmed_rows, "trimmed")
    if set(raw_by_key) != set(trimmed_by_key):
        raise ValueError("Raw and trimmed QC sample sets disagree")
    output: list[dict[str, str]] = []
    for key in sorted(raw_by_key):
        raw = raw_by_key[key]
        trimmed = trimmed_by_key[key]
        raw_count = int(raw["total_sequences"])
        trimmed_count = int(trimmed["total_sequences"])
        if raw_count <= 0 or not 0 <= trimmed_count <= raw_count:
            raise ValueError(f"Invalid raw/trimmed read counts for {key}")
        identity_fields = ("study_id", "bioproject", "sample_id")
        if any(raw[field] != trimmed[field] for field in identity_fields):
            raise ValueError(f"Raw/trimmed sample metadata disagree for {key}")
        output.append(
            {
                "study_id": raw["study_id"],
                "bioproject": raw["bioproject"],
                "sample_id": raw["sample_id"],
                "run_accession": raw["run_accession"],
                "read_direction": raw["read_direction"],
                "raw_total_sequences": str(raw_count),
                "trimmed_total_sequences": str(trimmed_count),
                "count_retention_percent": f"{trimmed_count / raw_count * 100:.6f}",
                "raw_sequence_length": raw["sequence_length"],
                "trimmed_sequence_length": trimmed["sequence_length"],
                "raw_per_base_quality_status": raw["per_base_quality_status"],
                "trimmed_per_base_quality_status": trimmed["per_base_quality_status"],
                "raw_adapter_content_status": raw["adapter_content_status"],
                "trimmed_adapter_content_status": trimmed["adapter_content_status"],
            }
        )
    return output


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
