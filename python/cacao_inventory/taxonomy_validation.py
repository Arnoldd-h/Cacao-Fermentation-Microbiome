"""Independent checks of R taxonomy masks, coverage and ASV conservation."""

from __future__ import annotations

import re
from pathlib import Path

from .io import read_tsv


def keyed(rows: list[dict], columns: tuple[str, ...]) -> dict[tuple[str, ...], dict]:
    result = {tuple(row[key] for key in columns): row for row in rows}
    if len(result) != len(rows) or any(any(not key for key in keys) for keys in result):
        raise ValueError(f"Missing or duplicate keys: {columns}")
    return result


def nonnegative_integer(value: str) -> int:
    if not re.fullmatch(r"[0-9]+", value):
        raise ValueError(f"Invalid integer: {value!r}")
    return int(value)


def mask_lineage(raw: dict, boots: dict, ranks: list[str], threshold: int) -> dict[str, str]:
    supported = True
    result = {}
    for rank in ranks:
        boot = nonnegative_integer(boots[rank])
        if boot > 100:
            raise ValueError("Bootstrap support outside 0..100")
        supported = supported and bool(raw[rank]) and boot >= threshold
        result[rank] = raw[rank] if supported else ""
    return result


def screening_flags(lineage: dict, settings: dict) -> dict[str, str]:
    kingdom = lineage["Kingdom"]
    status = "unclassified" if not kingdom else "bacterial" if kingdom == settings["bacterial_kingdom"] else "non_bacterial"
    labels = {value.casefold() for value in lineage.values() if value}
    organs = ";".join(value for value in settings["organelle_labels"] if value.casefold() in labels)
    named = bool(lineage["Genus"]) and not re.search(settings["unresolved_genus_pattern"], lineage["Genus"], re.I)
    return {"kingdom_status": status, "organelle_flag": organs,
            "named_genus": "TRUE" if named else "FALSE", "excluded": "FALSE"}


def validate_taxonomy_tables(directory: Path, input_directory: Path, config: dict) -> dict:
    ranks = config["classification"]["tax_levels"]
    primary = config["classification"]["primary_min_boot"]
    thresholds = [primary, *config["classification"]["sensitivity_min_boot"]]
    identifiers = ["study_id", "asv_id"]
    source_rows = read_tsv(input_directory / "asv_sequences.tsv", identifiers + ["sequence", "total_reads"])
    sequences = keyed(source_rows, ("study_id", "asv_id"))
    studies = {key[0] for key in sequences}
    if len(studies) != 1 or not sequences:
        raise ValueError("Taxonomy requires exactly one nonempty study")
    study = next(iter(studies))
    counts_rows = read_tsv(input_directory / "asv_counts.tsv", ["study_id", "sample_id"])
    counts = keyed(counts_rows, ("study_id", "sample_id"))
    if not counts or any(key[0] != study for key in counts):
        raise ValueError("Sample counts mix studies")
    totals = {key: 0 for key in sequences}
    for row in counts.values():
        if set(row) - {"study_id", "sample_id"} != {key[1] for key in sequences}:
            raise ValueError("Count columns disagree with ASVs")
        for key in sequences:
            totals[key] += nonnegative_integer(row[key[1]])
    for key, row in sequences.items():
        if totals[key] != nonnegative_integer(row["total_reads"]) or not totals[key]:
            raise ValueError("ASV totals disagree with original counts")
    raw = keyed(read_tsv(directory / "taxonomy_unfiltered.tsv", identifiers + ranks), tuple(identifiers))
    boots = keyed(read_tsv(directory / "taxonomy_bootstraps.tsv", identifiers + ranks), tuple(identifiers))
    assigned = keyed(read_tsv(directory / "taxonomy.tsv", identifiers + ranks), tuple(identifiers))
    flags = keyed(read_tsv(directory / "taxonomy_screening.tsv", identifiers + ["total_reads", "excluded"]), tuple(identifiers))
    if any(set(table) != set(sequences) for table in (raw, boots, assigned, flags)):
        raise ValueError("Taxonomy membership differs from original ASVs")
    sensitivity = keyed(read_tsv(directory / "taxonomy_sensitivity.tsv", identifiers + ["min_boot"] + ranks),
                        ("study_id", "asv_id", "min_boot"))
    expected_sensitivity = {(study, key[1], str(t)) for key in sequences for t in thresholds}
    if set(sensitivity) != expected_sensitivity:
        raise ValueError("Sensitivity ASV/threshold membership mismatch")
    expected = {}
    for key in sequences:
        for threshold in thresholds:
            lineage = mask_lineage(raw[key], boots[key], ranks, threshold)
            expected[key, threshold] = lineage
            if any(sensitivity[(*key, str(threshold))][rank] != lineage[rank] for rank in ranks):
                raise ValueError("Sensitivity taxonomy disagrees with bootstrap mask")
        primary_lineage = expected[key, primary]
        if any(assigned[key][rank] != primary_lineage[rank] for rank in ranks):
            raise ValueError("Primary taxonomy disagrees with bootstrap mask")
        computed_flags = screening_flags(primary_lineage, config["screening"])
        if any(flags[key].get(field) != value for field, value in computed_flags.items()):
            raise ValueError("Screening flags disagree with taxonomy or exclusions occurred")
        if nonnegative_integer(flags[key]["total_reads"]) != totals[key]:
            raise ValueError("Screening changed read totals")
    coverage = keyed(read_tsv(directory / "assignment_coverage.tsv", ["study_id", "min_boot", "rank", "assigned_asvs", "total_asvs", "assigned_reads", "total_reads"]),
                     ("study_id", "min_boot", "rank"))
    expected_keys = {(study, str(t), rank) for t in thresholds for rank in ranks}
    if set(coverage) != expected_keys:
        raise ValueError("Coverage rank/threshold membership mismatch")
    for threshold in thresholds:
        for rank in ranks:
            selected = [key for key in sequences if expected[key, threshold][rank]]
            row = coverage[study, str(threshold), rank]
            calculated = {"assigned_asvs": len(selected), "total_asvs": len(sequences),
                          "assigned_reads": sum(totals[key] for key in selected), "total_reads": sum(totals.values())}
            if any(nonnegative_integer(row[field]) != value for field, value in calculated.items()):
                raise ValueError("Coverage totals disagree with original reads")
    sample_rows = keyed(read_tsv(directory / "sample_coverage.tsv", ["study_id", "sample_id", "min_boot", "rank", "assigned_reads", "total_reads"]),
                        ("study_id", "sample_id", "min_boot", "rank"))
    expected_keys = {(study, key[1], str(t), rank) for key in counts for t in thresholds for rank in ranks}
    if set(sample_rows) != expected_keys:
        raise ValueError("Sample coverage membership mismatch")
    for key, count_row in counts.items():
        total = sum(nonnegative_integer(count_row[asv[1]]) for asv in sequences)
        for threshold in thresholds:
            for rank in ranks:
                selected = [asv for asv in sequences if expected[asv, threshold][rank]]
                retained = sum(nonnegative_integer(count_row[asv[1]]) for asv in selected)
                row = sample_rows[(*key, str(threshold), rank)]
                if nonnegative_integer(row["assigned_reads"]) != retained or nonnegative_integer(row["total_reads"]) != total:
                    raise ValueError("Sample coverage disagrees with original reads")
    return {"study_id": study, "asvs": len(sequences), "samples": len(counts), "total_reads": sum(totals.values()),
            "primary_min_boot": primary, "sensitivity_min_boot": thresholds[1:],
            "genus_assigned_asvs": sum(bool(expected[key, primary]["Genus"]) for key in sequences),
            "named_genus_asvs": sum(flags[key]["named_genus"] == "TRUE" for key in sequences),
            "organelle_flagged_asvs": sum(bool(flags[key]["organelle_flag"]) for key in sequences),
            "non_bacterial_asvs": sum(flags[key]["kingdom_status"] == "non_bacterial" for key in sequences),
            "unclassified_kingdom_asvs": sum(flags[key]["kingdom_status"] == "unclassified" for key in sequences),
            "excluded_asvs": 0}
