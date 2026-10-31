"""Synthetic conservation/lineage cases; these are not biological observations."""

import copy
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
SPEC = importlib.util.spec_from_file_location("summarize_composition", ROOT / "scripts/reporting/summarize_composition.py")
REPORTING = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORTING)
from cacao_inventory.io import write_tsv_atomic


class CompositionReportingTests(unittest.TestCase):
    def setUp(self):
        self.counts = [{"study_id": "study", "sample_id": "sample", "a": "20", "b": "30", "c": "50", "d": "0"}]
        self.features = ["a", "b", "c", "d"]
        self.metadata = [{field: "" for field in REPORTING.IDENTITY}]
        self.metadata[0].update(study_id="study", sample_id="sample", run_accession="SRR1", fermentation_batch="batch", fermentation_hours="0", relative_time="0", fermentation_stage="early", sampling_stratum="1")
        self.taxonomy = []
        for asv, family, genus in (("a", "F1", "shared"), ("b", "F2", "shared"), ("c", "F1", ""), ("d", "F1", "absent")):
            self.taxonomy.append(dict(study_id="study", asv_id=asv, **dict(zip(REPORTING.LEVELS, ["Bacteria", "P", "C", "O", family, genus]))))
        self.ranks = ["Family", "Genus"]

    def summarize(self):
        return REPORTING.aggregate_taxa(self.counts, self.features, self.taxonomy, self.metadata, self.ranks)

    def test_homonyms_preserve_ancestors_and_unknowns_keep_full_denominator(self):
        composition, coverage = self.summarize()
        families = [row for row in composition if row["taxonomic_rank"] == "Family"]
        genera = [row for row in composition if row["taxonomic_rank"] == "Genus"]
        self.assertEqual({row["Family"]: int(row["read_count"]) for row in families}, {"F1": 70, "F2": 30})
        self.assertEqual({(row["Family"], row["Genus"]): int(row["read_count"]) for row in genera}, {("F1", ""): 50, ("F1", "shared"): 20, ("F2", "shared"): 30})
        self.assertEqual(sum(float(row["relative_abundance"]) for row in genera), 1)
        self.assertEqual(coverage[1]["labeled_reads"], "50")
        self.assertEqual(coverage[1]["unlabeled_reads"], "50")
        self.assertEqual(float(coverage[1]["labeled_fraction"]), 0.5)
        self.assertNotIn("absent", {row["Genus"] for row in genera})

    def test_duplicate_taxonomy_and_missing_metadata_fail(self):
        self.taxonomy.append(copy.deepcopy(self.taxonomy[0]))
        with self.assertRaises(ValueError):
            self.summarize()
        self.taxonomy.pop()
        self.metadata.clear()
        with self.assertRaises(ValueError):
            self.summarize()

    def test_mixed_studies_empty_samples_and_fractional_counts_fail(self):
        original = copy.deepcopy(self.counts)
        self.counts.append({**self.counts[0], "study_id": "other"})
        with self.assertRaises(ValueError):
            self.summarize()
        self.counts = [{**original[0], **{feature: "0" for feature in self.features}}]
        with self.assertRaises(ValueError):
            self.summarize()
        self.counts = [{**original[0], "a": "1.5"}]
        with self.assertRaises(ValueError):
            self.summarize()

    def test_source_order_does_not_change_taxonomic_rows(self):
        expected = self.summarize()
        self.features.reverse()
        self.taxonomy.reverse()
        self.assertEqual(self.summarize(), expected)

    def test_unclassified_kingdom_is_retained_without_invented_labels(self):
        self.taxonomy[2].update({level: "" for level in REPORTING.LEVELS})
        composition, coverage = self.summarize()
        unknown = [row for row in composition if row["taxonomic_rank"] == "Genus" and not row["Kingdom"]]
        self.assertEqual(len(unknown), 1)
        self.assertEqual(unknown[0]["read_count"], "50")
        self.assertTrue(all(unknown[0][level] == "" for level in REPORTING.LEVELS))
        self.assertEqual(coverage[1]["unlabeled_reads"], "50")

    def test_policies_reject_imputation_pooling_and_unsafe_outputs(self):
        config = json.loads((ROOT / "config/taxonomic_reporting.yaml").read_text())
        REPORTING.reporting_profile(config, "full_study", ROOT)
        for field, value in (("inference", True), ("sample_aggregation", "batch"), ("preserve_ancestral_ranks", False), ("unclassified_policy", "drop"), ("ranks", ["Genus", "Genus"])):
            altered = copy.deepcopy(config)
            altered[field] = value
            with self.assertRaises(ValueError):
                REPORTING.reporting_profile(altered, "full_study", ROOT)
        config["analyses"]["full_study"]["output_dir"] = "data/raw/output"
        with self.assertRaises(ValueError):
            REPORTING.reporting_profile(config, "full_study", ROOT)

    def test_validator_detects_altered_counts_metadata_denominator_and_extra_rows(self):
        self.counts[0]["d"] = "5"
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source"
            report = directory / "report"
            report.mkdir()
            write_tsv_atomic(source / "asv_counts.tsv", self.counts, ["study_id", "sample_id", *self.features])
            write_tsv_atomic(source / "taxonomy.tsv", self.taxonomy, ["study_id", "asv_id", *REPORTING.LEVELS])
            write_tsv_atomic(source / "sample_metadata.tsv", self.metadata, REPORTING.IDENTITY)
            composition, coverage = self.summarize()
            summary = {"study_id": "study", "samples": 1, "asvs": 4, "retained_reads": 105, "taxonomic_ranks": self.ranks, "composition_rows": 6, "sample_aggregation": "none", "inference": False}
            (report / "summary.json").write_text(json.dumps(summary))
            def write_outputs(rows):
                write_tsv_atomic(report / REPORTING.PRODUCTS[0], rows, REPORTING.COMPOSITION_COLUMNS)
                write_tsv_atomic(report / REPORTING.PRODUCTS[1], coverage, REPORTING.COVERAGE_COLUMNS)
            write_outputs(composition)
            paths = {"bacterial_dir": source}
            self.assertEqual(REPORTING.validate_tables(paths, report, self.ranks, 1, 1e-8), summary)
            for field, value in (("read_count", "71"), ("relative_abundance", "1.0"), ("fermentation_hours", "12"), ("Family", "wrong")):
                altered = copy.deepcopy(composition)
                altered[0][field] = value
                write_outputs(altered)
                with self.assertRaises(ValueError):
                    REPORTING.validate_tables(paths, report, self.ranks, 1, 1e-8)
            write_outputs(composition + [copy.deepcopy(composition[0])])
            with self.assertRaises(ValueError):
                REPORTING.validate_tables(paths, report, self.ranks, 1, 1e-8)


if __name__ == "__main__":
    unittest.main()
