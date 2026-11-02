"""Streaming checks for terminal-only paired-read trimming and FastQC coverage.

Length thresholds are descriptive diagnostics; these functions never filter reads.
"""

from collections import Counter
from contextlib import ExitStack
import gzip
from itertools import zip_longest
import math
from pathlib import Path
import re
from zipfile import ZipFile
import zlib


QUALITY_COLUMNS = (
    "Base", "Mean", "Median", "Lower Quartile", "Upper Quartile",
    "10th Percentile", "90th Percentile",
)
_NON_PRINTABLE = re.compile(rb"[^!-~]")
_BASE_BIN = re.compile(r"([1-9][0-9]*)(?:-([1-9][0-9]*))?\Z")


def _pair_id(header: bytes, direction: int) -> bytes:
    identifier = header[1:].split()[0]
    if identifier.endswith((b"/1", b"/2")):
        if not identifier.endswith(f"/{direction}".encode("ascii")):
            raise ValueError(f"FASTQ mate suffix disagrees with R{direction}")
        identifier = identifier[:-2]
    if not identifier:
        raise ValueError("FASTQ identifier must not be empty")
    return identifier


def _fastq_records(handle, path: Path):
    number = 0
    while True:
        header_line = handle.readline()
        if not header_line:
            return
        number += 1
        lines = [header_line, handle.readline(), handle.readline(), handle.readline()]
        if any(not line for line in lines):
            raise ValueError(f"Truncated FASTQ record {number}: {path}")
        header, sequence, separator, quality = [line.rstrip(b"\r\n") for line in lines]
        if not header.startswith(b"@") or not header[1:].split():
            raise ValueError(f"Invalid FASTQ header at record {number}: {path}")
        if not separator.startswith(b"+"):
            raise ValueError(f"Invalid FASTQ separator at record {number}: {path}")
        if len(sequence) != len(quality):
            raise ValueError(f"FASTQ sequence/quality length mismatch at record {number}: {path}")
        if _NON_PRINTABLE.search(sequence) or _NON_PRINTABLE.search(quality):
            raise ValueError(f"Invalid FASTQ sequence/quality character at record {number}: {path}")
        yield header, sequence, separator, quality


def _direction_counts():
    return {
        "changed_reads": 0, "bases_removed": 0, "below_min_len": 0,
        "baseline_below_min_len": 0, "newly_below_min_len": 0,
        "length_counts": Counter(), "baseline_length_counts": Counter(),
    }


def _audit_read(before, after, primer_length: int, counts: dict, min_len: int):
    if before[0] != after[0] or before[2] != after[2]:
        raise ValueError("FASTQ header or separator changed during terminal trimming")
    old_length, new_length = len(before[1]), len(after[1])
    removed = old_length - new_length
    if removed not in (0, primer_length):
        raise ValueError(f"Unexpected terminal trimming length: {removed}, expected 0 or {primer_length}")
    if before[1][:new_length] != after[1] or before[3][:new_length] != after[3]:
        raise ValueError("Terminal trimming altered retained sequence or quality")
    counts["changed_reads"] += removed > 0
    counts["bases_removed"] += removed
    counts["below_min_len"] += new_length < min_len
    counts["baseline_below_min_len"] += old_length < min_len
    counts["newly_below_min_len"] += old_length >= min_len > new_length
    counts["length_counts"][new_length] += 1
    counts["baseline_length_counts"][old_length] += 1


def _finish_counts(counts: dict):
    for field, prefix in (("length_counts", ""), ("baseline_length_counts", "baseline_")):
        histogram = counts[field]
        counts[f"{prefix}min_length"] = min(histogram) if histogram else None
        counts[f"{prefix}max_length"] = max(histogram) if histogram else None
        counts[field] = {str(length): histogram[length] for length in sorted(histogram)}


def audit_pair(inputs: list[Path], outputs: list[Path], primer_lengths: list[int], min_len: int) -> dict:
    """Verify four gzip FASTQs in lockstep and return JSON-compatible diagnostics.

    Each output must retain every input identifier, in order, and exactly the
    sequence/quality prefix obtained by removing zero or one full primer length.
    Empty sequences remain valid records. Corrupt gzip, malformed FASTQ, changed
    bases/qualities, altered names, mate desynchronization and lost reads fail.
    """
    if len(inputs) != 2 or len(outputs) != 2 or len(primer_lengths) != 2:
        raise ValueError("Exactly two inputs, outputs and primer lengths are required")
    if any(type(length) is not int or length <= 0 for length in primer_lengths):
        raise ValueError("Primer lengths must be positive integers")
    if type(min_len) is not int or min_len < 0:
        raise ValueError("The diagnostic minimum length must be a nonnegative integer")
    result = {
        "pairs": 0, "pairs_below_min_len": 0, "pairs_below_min_len_baseline": 0,
        "pairs_newly_below_min_len": 0,
        "directions": {"R1": _direction_counts(), "R2": _direction_counts()},
    }
    missing = object()
    try:
        with ExitStack() as stack:
            paths = [Path(path) for path in [*inputs, *outputs]]
            readers = [_fastq_records(stack.enter_context(gzip.open(path, "rb")), path) for path in paths]
            for records in zip_longest(*readers, fillvalue=missing):
                if any(record is missing for record in records):
                    raise ValueError("FASTQ input/output mates contain different numbers of reads")
                old_r1, old_r2, new_r1, new_r2 = records
                if _pair_id(old_r1[0], 1) != _pair_id(old_r2[0], 2):
                    raise ValueError("FASTQ mates have different identifiers or order")
                for index, (before, after) in enumerate(((old_r1, new_r1), (old_r2, new_r2))):
                    _audit_read(before, after, primer_lengths[index], result["directions"][f"R{index + 1}"], min_len)
                old_below = min(len(old_r1[1]), len(old_r2[1])) < min_len
                new_below = min(len(new_r1[1]), len(new_r2[1])) < min_len
                result["pairs"] += 1
                result["pairs_below_min_len_baseline"] += old_below
                result["pairs_below_min_len"] += new_below
                result["pairs_newly_below_min_len"] += new_below and not old_below
    except (OSError, EOFError, zlib.error) as exc:
        raise ValueError(f"Cannot read complete gzip FASTQ inputs/outputs: {exc}") from exc
    for counts in result["directions"].values():
        _finish_counts(counts)
    return result


def _parse_modules(lines: list[str]) -> tuple[dict, dict]:
    if not lines or not re.fullmatch(r"##FastQC\t\S+", lines[0]):
        raise ValueError("FastQC version header is missing or invalid")
    modules, status = {}, {}
    current = None
    for line in lines[1:]:
        if line == ">>END_MODULE":
            if current is None:
                raise ValueError("FastQC END_MODULE outside a module")
            current = None
        elif line.startswith(">>"):
            fields = line[2:].split("\t")
            if current is not None or len(fields) != 2 or fields[1] not in ("pass", "warn", "fail"):
                raise ValueError("Incomplete or invalid FastQC module")
            current, module_status = fields
            if not current or current in modules:
                raise ValueError("Empty or duplicated FastQC module")
            modules[current], status[current] = [], module_status
        elif current is None:
            if line:
                raise ValueError("FastQC data outside a module")
        else:
            modules[current].append(line)
    if current is not None:
        raise ValueError("FastQC module is missing END_MODULE")
    for module in ("Basic Statistics", "Per base sequence quality"):
        if module not in modules:
            raise ValueError(f"Required FastQC module is missing: {module}")
    return modules, status


def _parse_base_bin(value: str) -> tuple[int, int]:
    match = _BASE_BIN.fullmatch(value)
    if match is None:
        raise ValueError(f"Invalid FastQC base bin: {value!r}")
    start = int(match.group(1))
    end = int(match.group(2) or start)
    if end < start:
        raise ValueError(f"Reversed FastQC base bin: {value!r}")
    return start, end


def _quality_rows(lines: list[str]) -> list[dict]:
    if not lines or lines[0] != "#" + "\t".join(QUALITY_COLUMNS):
        raise ValueError("Unexpected FastQC per-base quality columns")
    rows = []
    previous_end = 0
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(QUALITY_COLUMNS):
            raise ValueError("Malformed FastQC per-base quality row")
        row = dict(zip(QUALITY_COLUMNS, fields))
        start, end = _parse_base_bin(row["Base"])
        if start != previous_end + 1:
            raise ValueError("FastQC quality bins must be contiguous from base 1")
        previous_end = end
        for key in QUALITY_COLUMNS[1:]:
            # FastQC reports unavailable percentiles as literal NaN in sparse
            # bins; retain that text rather than inventing a zero or dropping it.
            if row[key] != "NaN" and not math.isfinite(float(row[key])):
                raise ValueError(f"Invalid nonfinite FastQC quality: {key}")
        rows.append(row)
    return rows


def read_fastqc(path: Path) -> dict:
    """Read one FastQC ZIP, requiring complete modules and standard QC headers.

    Basic statistics and quality values retain FastQC's original text columns.
    All module statuses are retained, including warnings and failures. Native
    "NaN" quality values remain text denoting unavailable statistics.
    """
    with ZipFile(path) as archive:
        members = [entry for entry in archive.infolist() if entry.filename.endswith("/fastqc_data.txt")]
        if len(members) != 1:
            raise ValueError("FastQC ZIP must contain exactly one /fastqc_data.txt")
        lines = archive.read(members[0]).decode("utf-8").splitlines()
    modules, status = _parse_modules(lines)
    basic_lines = modules["Basic Statistics"]
    if not basic_lines or basic_lines[0] != "#Measure\tValue":
        raise ValueError("Unexpected FastQC basic-statistics columns")
    basic = {}
    for line in basic_lines[1:]:
        fields = line.split("\t")
        if len(fields) != 2 or not fields[0] or fields[0] in basic:
            raise ValueError("Malformed or duplicated FastQC basic statistic")
        basic[fields[0]] = fields[1]
    total = basic.get("Total Sequences", "")
    if not re.fullmatch(r"[0-9]+", total):
        raise ValueError("FastQC Total Sequences must be a nonnegative integer")
    rows = _quality_rows(modules["Per base sequence quality"])
    return {"basic": basic, "module_status": status, "quality_rows": rows}


def quality_coverage(quality_rows: list[dict], length_counts: dict) -> list[dict]:
    """Annotate 1-based inclusive FastQC bins with exact read/base coverage.

    The denominator includes every read, including zero-length reads. A bin's
    base observations sum the number of retained bases in that bin over reads;
    this is not a reconstruction or reweighting of FastQC's quality statistic.
    Bins must cover every observed position continuously from base 1. The final
    bin may extend past the longest read, but no whole bin may lack coverage.
    An empty quality report is valid only when the histogram contains no bases.
    """
    histogram = {}
    for key, count in length_counts.items():
        if isinstance(key, bool) or not re.fullmatch(r"0|[1-9][0-9]*", str(key)):
            raise ValueError("Read lengths must be nonnegative integer histogram keys")
        if type(count) is not int or count < 0 or int(key) in histogram:
            raise ValueError("Read length counts must be unique nonnegative integers")
        histogram[int(key)] = count
    total = sum(histogram.values())
    max_length = max((length for length, count in histogram.items() if count > 0), default=0)
    rows = []
    previous_end = 0
    for row in quality_rows:
        start, end = _parse_base_bin(row["Base"])
        if start != previous_end + 1:
            raise ValueError("FastQC quality bins must be contiguous from base 1")
        if start > max_length:
            raise ValueError("FastQC quality bin contains no observed bases")
        previous_end = end
        reaching_start = sum(count for length, count in histogram.items() if length >= start)
        reaching_end = sum(count for length, count in histogram.items() if length >= end)
        observations = sum(max(0, min(length, end) - start + 1) * count for length, count in histogram.items())
        rows.append({
            **row, "bin_start": start, "bin_end": end,
            "reads_reaching_start": reaching_start, "reads_reaching_end": reaching_end,
            "base_observations": observations,
            "fraction_reads_reaching_start": reaching_start / total if total else None,
        })
    if previous_end < max_length:
        raise ValueError("FastQC quality report omits observed terminal bases")
    return rows
