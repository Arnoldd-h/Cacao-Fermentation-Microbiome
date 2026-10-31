#!/usr/bin/env python3
"""Summarize validated within-sample counts by their original taxonomic lineage.

Inputs: reporting policy/profile and independently validated bacterial tables.
Outputs: sample_composition.tsv, sample_label_coverage.tsv, summary.json, hashes,
SUCCESS and validation.json under the registered composition directory.
Fails on unsupported policies, mixed studies, missing/duplicate identities,
altered upstream data, or counts/labels/metadata that do not reconcile.
--validate-only checks existing outputs without regenerating them.
No sample aggregation, taxonomy imputation, pseudocounts or inferential tests.
"""

import argparse
import json
import logging
import math
import platform
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/filtering"))
from validate_bacterial_table import validate as validate_filtering
from cacao_inventory.analysis_artifacts import begin_stage, finish_stage, read_counts, read_table, unique_rows, validate_artifacts, validate_policy, write_json
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.provenance import file_records
from cacao_inventory.taxonomy_reference import repository_path

LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
IDENTITY = ["study_id", "sample_id", "run_accession", "fermentation_batch", "fermentation_hours", "relative_time", "fermentation_stage", "sampling_stratum"]
COMPOSITION_COLUMNS = [*IDENTITY, "taxonomic_rank", *LEVELS, "read_count", "total_retained_reads", "relative_abundance"]
COVERAGE_COLUMNS = [*IDENTITY, "taxonomic_rank", "labeled_reads", "unlabeled_reads", "total_retained_reads", "labeled_fraction"]
PRODUCTS = ["sample_composition.tsv", "sample_label_coverage.tsv", "summary.json"]


def reporting_profile(config, analysis, root):
    expected = {"preserve_ancestral_ranks": True, "unclassified_policy": "preserve_empty_lineage",
                "sample_aggregation": "none", "zero_policy": "omit_zero_count_rows",
                "relative_abundance_denominator": "all_retained_reads_per_sample", "inference": False}
    if any(config.get(key) != value or (isinstance(value, bool) and type(config.get(key)) is not bool) for key, value in expected.items()):
        raise ValueError("Unsupported taxonomic reporting policy")
    ranks = config.get("ranks", [])
    if not ranks or len(set(ranks)) != len(ranks) or any(rank not in ("Family", "Genus") for rank in ranks):
        raise ValueError("Reporting requires unique Family/Genus ranks")
    profile = config["analyses"].get(analysis)
    if not profile or type(profile["expected_samples"]) is not int or profile["expected_samples"] < 1:
        raise ValueError("Unknown reporting analysis or invalid sample count")
    output = repository_path(root, profile["output_dir"])
    if not output.is_relative_to(root.resolve() / "results/composition"):
        raise ValueError("Composition outputs must remain under results/composition")
    source = repository_path(root, profile["diversity_config"])
    return profile, source, output


def load_sources(paths):
    counts, features = read_counts(paths["bacterial_dir"] / "asv_counts.tsv")
    taxonomy = read_table(paths["bacterial_dir"] / "taxonomy.tsv", ["study_id", "asv_id", *LEVELS])
    metadata = read_table(paths["bacterial_dir"] / "sample_metadata.tsv", IDENTITY)
    return counts, features, taxonomy, metadata


def aggregate_taxa(counts, features, taxonomy, metadata, ranks):
    samples = unique_rows(counts, ("study_id", "sample_id"))
    calls = unique_rows(taxonomy, ("study_id", "asv_id"))
    meta = unique_rows(metadata, ("study_id", "sample_id"))
    studies = {row["study_id"] for row in counts}
    if len(studies) != 1 or not features or len(set(features)) != len(features):
        raise ValueError("Reporting requires one nonempty study and unique ASVs")
    study = next(iter(studies))
    if set(meta) != set(samples) or set(calls) != {(study, feature) for feature in features}:
        raise ValueError("Counts, taxonomy and metadata memberships differ")
    composition, coverage = [], []
    for key, sample in samples.items():
        if any(not isinstance(sample[feature], str) or not sample[feature].isdigit() for feature in features):
            raise ValueError("Counts must be nonnegative integer strings")
        total = sum(int(sample[feature]) for feature in features)
        if total == 0:
            raise ValueError("Empty sample in retained counts")
        identity = {field: meta[key][field] for field in IDENTITY}
        for rank in ranks:
            ancestors = LEVELS[:LEVELS.index(rank) + 1]
            groups = Counter()
            for feature in features:
                number = int(sample[feature])
                if number:
                    groups[tuple(calls[study, feature][level] for level in ancestors)] += number
            for lineage, number in sorted(groups.items()):
                labels = {level: "" for level in LEVELS}
                labels.update(zip(ancestors, lineage))
                composition.append({**identity, "taxonomic_rank": rank, **labels, "read_count": str(number),
                                    "total_retained_reads": str(total), "relative_abundance": repr(number / total)})
            labeled = sum(number for lineage, number in groups.items() if lineage[-1])
            coverage.append({**identity, "taxonomic_rank": rank, "labeled_reads": str(labeled),
                             "unlabeled_reads": str(total - labeled), "total_retained_reads": str(total),
                             "labeled_fraction": repr(labeled / total)})
    return composition, coverage


def validate_tables(paths, report, ranks, expected_samples, tolerance):
    counts, features, calls, metadata = load_sources(paths)
    samples = unique_rows(counts, ("study_id", "sample_id"))
    taxonomy = unique_rows(calls, ("study_id", "asv_id"))
    meta = unique_rows(metadata, ("study_id", "sample_id"))
    if len(samples) != expected_samples or set(meta) != set(samples):
        raise ValueError("Reporting sample scope changed")
    observed = read_table(report / PRODUCTS[0], COMPOSITION_COLUMNS)
    observed_keys = ("study_id", "sample_id", "taxonomic_rank", *LEVELS)
    rows = {}
    for row in observed:
        key = tuple(row[field] for field in observed_keys)
        if key in rows or any(not row[field] for field in ("study_id", "sample_id", "taxonomic_rank")):
            raise ValueError("Empty sample identity or duplicate taxonomic lineage")
        rows[key] = row
    observed_coverage = unique_rows(read_table(report / PRODUCTS[1], COVERAGE_COLUMNS), ("study_id", "sample_id", "taxonomic_rank"))
    expected_membership = set()
    for key, sample in samples.items():
        total = sum(int(sample[feature]) for feature in features)
        for rank in ranks:
            end = LEVELS.index(rank) + 1
            expected = {}
            for feature in features:
                if int(sample[feature]) == 0:
                    continue
                lineage = tuple(taxonomy[key[0], feature][level] if index < end else "" for index, level in enumerate(LEVELS))
                row_key = (*key, rank, *lineage)
                expected[row_key] = expected.get(row_key, 0) + int(sample[feature])
            expected_membership.update(expected)
            labeled = 0
            for row_key, number in expected.items():
                row = rows.get(row_key)
                if row is None or any(row[field] != meta[key][field] for field in IDENTITY) or row["read_count"] != str(number) or row["total_retained_reads"] != str(total):
                    raise ValueError("Composition count, lineage or metadata mismatch")
                value = float(row["relative_abundance"])
                if not math.isfinite(value) or not math.isclose(value, number / total, rel_tol=0, abs_tol=tolerance):
                    raise ValueError("Composition denominator or proportion mismatch")
                if row[rank]:
                    labeled += number
            cover = observed_coverage.get((*key, rank))
            if cover is None or any(cover[field] != meta[key][field] for field in IDENTITY) or any(cover[field] != str(value) for field, value in (("labeled_reads", labeled), ("unlabeled_reads", total - labeled), ("total_retained_reads", total))):
                raise ValueError("Label coverage count or metadata mismatch")
            fraction = float(cover["labeled_fraction"])
            if not math.isfinite(fraction) or not math.isclose(fraction, labeled / total, rel_tol=0, abs_tol=tolerance):
                raise ValueError("Label coverage fraction mismatch")
    if set(rows) != expected_membership or set(observed_coverage) != {(*key, rank) for key in samples for rank in ranks}:
        raise ValueError("Unexpected or missing composition/coverage rows")
    summary = {"study_id": counts[0]["study_id"], "samples": len(samples), "asvs": len(features),
               "retained_reads": sum(int(row[feature]) for row in counts for feature in features),
               "taxonomic_ranks": ranks, "composition_rows": len(observed), "sample_aggregation": "none", "inference": False}
    if json.loads((report / PRODUCTS[2]).read_text()) != summary:
        raise ValueError("Composition summary mismatch")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/taxonomic_reporting.yaml")
    parser.add_argument("--analysis", required=True, help="Analysis key in the reporting configuration")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    config = load_json_yaml(args.config)
    profile, source_config, report = reporting_profile(config, args.analysis, ROOT)
    diversity = load_json_yaml(source_config)
    paths = validate_policy(diversity, ROOT)
    upstream = validate_filtering(source_config, ROOT)
    if upstream["samples"] != profile["expected_samples"]:
        raise ValueError("Upstream sample scope differs from the reporting profile")
    tolerance = diversity["analysis"]["validation_tolerance"]
    if not args.validate_only:
        state = begin_stage(report, ROOT)
        inputs = [args.config, source_config, Path(__file__), ROOT / "scripts/filtering/validate_bacterial_table.py",
                  *[paths["bacterial_dir"] / name for name in ("asv_counts.tsv", "taxonomy.tsv", "sample_metadata.tsv")],
                  *[paths["filtering_dir"] / name for name in ("validation.json", "provenance.json")],
                  *[ROOT / "python/cacao_inventory" / (name + ".py") for name in ("analysis_artifacts", "bacterial_filter", "config", "io", "provenance", "taxonomy_reference")]]
        before = file_records(inputs, ROOT)
        source = load_sources(paths)
        composition, coverage = aggregate_taxa(*source, config["ranks"])
        write_tsv_atomic(report / PRODUCTS[0], composition, COMPOSITION_COLUMNS)
        write_tsv_atomic(report / PRODUCTS[1], coverage, COVERAGE_COLUMNS)
        counts, features = source[:2]
        write_json(report / PRODUCTS[2], {"study_id": counts[0]["study_id"], "samples": len(counts), "asvs": len(features),
                   "retained_reads": sum(int(row[feature]) for row in counts for feature in features),
                   "taxonomic_ranks": config["ranks"], "composition_rows": len(composition), "sample_aggregation": "none", "inference": False})
        finish_stage(report, ROOT, state, before, args.config, config, [report / name for name in PRODUCTS],
                     upstream, {"Python": platform.python_version()}, sys.argv)
    expected = [*[report / name for name in PRODUCTS], report / "config_snapshot.yaml", report / "input_checksums.json"]
    result = validate_artifacts(report, ROOT, expected, config)
    result.update(validate_tables(paths, report, config["ranks"], profile["expected_samples"], tolerance))
    if not args.validate_only:
        write_json(report / "validation.json", result)
    logging.info("Validated lineage-preserving composition for %d samples", result["samples"])
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
