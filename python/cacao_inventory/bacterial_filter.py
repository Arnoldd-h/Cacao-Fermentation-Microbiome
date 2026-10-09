"""Build auditable bacterial candidate counts without changing source tables."""

from pathlib import Path

from .analysis_artifacts import read_counts, read_table, unique_rows

DATA_PRODUCTS = ["asv_counts.tsv", "asv_sequences.tsv", "taxonomy.tsv", "sample_metadata.tsv"]
REPORT_PRODUCTS = ["asv_filter_log.tsv", "sample_retention.tsv", "summary.json"]


def exclusion_reason(lineage: dict, settings: dict) -> str:
    labels = {lineage[rank].casefold() for rank in settings["tax_levels"]}
    organs = [label for label in settings["organelle_labels"] if label.casefold() in labels]
    if organs:
        return "organelle:" + ";".join(organs)
    if settings["exclude_known_non_bacterial"] and lineage["Kingdom"] and lineage["Kingdom"] != settings["bacterial_kingdom"]:
        return "non_bacterial_kingdom:" + lineage["Kingdom"]
    return ""


def prepare_tables(dada2: Path, taxonomy: Path, settings: dict) -> tuple[dict, dict]:
    counts, features = read_counts(dada2 / "asv_counts.tsv")
    ranks = settings["tax_levels"]
    sequences = read_table(dada2 / "asv_sequences.tsv", ["study_id", "asv_id", "sequence", "total_reads"])
    calls = read_table(taxonomy / "taxonomy.tsv", ["study_id", "asv_id", *ranks])
    metadata = read_table(dada2 / "sample_metadata.tsv", ["study_id", "sample_id"])
    keyed_sequences = unique_rows(sequences, ("study_id", "asv_id"))
    keyed_calls = unique_rows(calls, ("study_id", "asv_id"))
    study = counts[0]["study_id"]
    if set(keyed_sequences) != {(study, feature) for feature in features} or set(keyed_calls) != set(keyed_sequences):
        raise ValueError("Counts, sequences and taxonomy membership disagree")
    if set(unique_rows(metadata, ("study_id", "sample_id"))) != {(study, row["sample_id"]) for row in counts}:
        raise ValueError("Sample metadata membership disagrees")
    reasons = {feature: exclusion_reason(keyed_calls[study, feature], settings) for feature in features}
    kept = [feature for feature in features if not reasons[feature]]
    if not kept:
        raise ValueError("Filtering leaves no ASVs")
    log = []
    for feature in features:
        total = sum(int(row[feature]) for row in counts)
        if str(total) != keyed_sequences[study, feature]["total_reads"]:
            raise ValueError("Sequence total differs from counts")
        call = keyed_calls[study, feature]
        log.append({"study_id": study, "asv_id": feature, "retained": "FALSE" if reasons[feature] else "TRUE",
            "reason": reasons[feature] or ("retained_bacterial" if call["Kingdom"] else "retained_unclassified_kingdom"),
            "total_reads": str(total), **{rank: call[rank] for rank in ranks}})
    retained_counts, retention = [], []
    for row in counts:
        total = sum(int(row[feature]) for feature in features)
        retained = sum(int(row[feature]) for feature in kept)
        if not retained:
            raise ValueError(f"Filtering empties sample {row['sample_id']}; no silent sample exclusion")
        retained_counts.append({key: row[key] for key in ("study_id", "sample_id", *kept)})
        retention.append({"study_id": study, "sample_id": row["sample_id"], "input_reads": str(total),
            "retained_reads": str(retained), "excluded_reads": str(total - retained),
            "retained_fraction": format(retained / total, ".12g")})
    summary = {"study_id": study, "samples": len(counts), "input_asvs": len(features), "retained_asvs": len(kept),
        "excluded_asvs": len(features) - len(kept), "input_reads": sum(int(row["input_reads"]) for row in retention),
        "retained_reads": sum(int(row["retained_reads"]) for row in retention),
        "excluded_reads": sum(int(row["excluded_reads"]) for row in retention),
        "unclassified_kingdom_retained_asvs": sum(not keyed_calls[study, feature]["Kingdom"] for feature in kept),
        "sample_exclusions": 0}
    return {"asv_counts.tsv": retained_counts, "asv_sequences.tsv": [row for row in sequences if row["asv_id"] in kept],
            "taxonomy.tsv": [row for row in calls if row["asv_id"] in kept], "sample_metadata.tsv": metadata,
            "asv_filter_log.tsv": log, "sample_retention.tsv": retention}, summary
