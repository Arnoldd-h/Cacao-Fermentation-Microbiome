"""Synthetic taxonomies exercise scientific validation; no biological evidence."""

import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.config import load_json_yaml
from cacao_inventory.taxonomy_validation import mask_lineage, screening_flags, validate_taxonomy_tables


def write_rows(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class TaxonomyTests(unittest.TestCase):
    def setUp(self):
        self.config = load_json_yaml(ROOT / "config/taxonomy.yaml")
        self.ranks = self.config["classification"]["tax_levels"]
        self.raw = dict(zip(self.ranks, ["Bacteria", "P", "C", "O", "F", "G"]))
        self.boot = dict(zip(self.ranks, ["100", "100", "100", "100", "100", "70"]))

    def test_threshold_masks_only_unsupported_descendants(self):
        self.assertEqual(mask_lineage(self.raw, self.boot, self.ranks, 80)["Genus"], "")
        self.assertEqual(mask_lineage(self.raw, self.boot, self.ranks, 50)["Genus"], "G")

    def test_low_support_or_missing_parent_removes_high_support_child(self):
        for raw, boot in (({**self.raw, "Family": ""}, self.boot), (self.raw, {**self.boot, "Family": "40", "Genus": "100"})):
            self.assertEqual(mask_lineage(raw, boot, self.ranks, 50)["Genus"], "")

    def test_bootstrap_outside_range_or_fractional_fails(self):
        for value in ("101", "-1", "50.1"):
            with self.assertRaises(ValueError):
                mask_lineage(self.raw, {**self.boot, "Genus": value}, self.ranks, 80)

    def test_organelle_uncultured_and_unknown_flags_preserve_asvs(self):
        flags = screening_flags({**self.raw, "Order": "Chloroplast", "Genus": "uncultured bacterium"}, self.config["screening"])
        self.assertEqual(flags["organelle_flag"], "Chloroplast")
        self.assertEqual(flags["named_genus"], "FALSE")
        self.assertEqual(flags["excluded"], "FALSE")
        flags = screening_flags(dict.fromkeys(self.ranks, ""), self.config["screening"])
        self.assertEqual(flags["kingdom_status"], "unclassified")

    def fixture(self, root):
        inputs, outputs = root / "inputs", root / "outputs"
        inputs.mkdir(); outputs.mkdir()
        identifiers = {"study_id": "fixture", "asv_id": "fixture__ASV1"}
        write_rows(inputs / "asv_sequences.tsv", [{**identifiers, "sequence": "ACGT", "total_reads": "6"}])
        write_rows(inputs / "asv_counts.tsv", [{"study_id": "fixture", "sample_id": "s1", "fixture__ASV1": "6"}])
        write_rows(outputs / "taxonomy_unfiltered.tsv", [{**identifiers, **self.raw}])
        write_rows(outputs / "taxonomy_bootstraps.tsv", [{**identifiers, **self.boot}])
        primary = mask_lineage(self.raw, self.boot, self.ranks, 80)
        write_rows(outputs / "taxonomy.tsv", [{**identifiers, **primary}])
        write_rows(outputs / "taxonomy_screening.tsv", [{**identifiers, "total_reads": "6", **screening_flags(primary, self.config["screening"])}])
        sensitivity, coverage, samples = [], [], []
        for threshold in (80, 50):
            masked = mask_lineage(self.raw, self.boot, self.ranks, threshold)
            sensitivity.append({**identifiers, "min_boot": str(threshold), **masked})
            for rank in self.ranks:
                assigned = bool(masked[rank])
                coverage.append({"study_id": "fixture", "min_boot": str(threshold), "rank": rank,
                                 "assigned_asvs": str(int(assigned)), "total_asvs": "1", "assigned_reads": "6" if assigned else "0", "total_reads": "6"})
                samples.append({"study_id": "fixture", "sample_id": "s1", "min_boot": str(threshold), "rank": rank,
                                "assigned_reads": "6" if assigned else "0", "total_reads": "6"})
        write_rows(outputs / "taxonomy_sensitivity.tsv", sensitivity)
        write_rows(outputs / "assignment_coverage.tsv", coverage)
        write_rows(outputs / "sample_coverage.tsv", samples)
        return inputs, outputs

    def test_full_table_conservation(self):
        with tempfile.TemporaryDirectory() as temporary:
            inputs, outputs = self.fixture(Path(temporary))
            summary = validate_taxonomy_tables(outputs, inputs, self.config)
            self.assertEqual(summary["total_reads"], 6)
            self.assertEqual(summary["genus_assigned_asvs"], 0)

    def test_unjustified_genus_call_fails_even_with_consistent_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            inputs, outputs = self.fixture(Path(temporary))
            write_rows(outputs / "taxonomy.tsv", [{"study_id": "fixture", "asv_id": "fixture__ASV1", **self.raw}])
            with self.assertRaisesRegex(ValueError, "bootstrap mask"):
                validate_taxonomy_tables(outputs, inputs, self.config)

    def test_lost_asv_or_exclusion_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            inputs, outputs = self.fixture(Path(temporary))
            path = outputs / "taxonomy_screening.tsv"
            path.write_text(path.read_text().replace("FALSE\tFALSE", "FALSE\tTRUE"))
            with self.assertRaisesRegex(ValueError, "exclusions"):
                validate_taxonomy_tables(outputs, inputs, self.config)

    def test_coverage_cannot_change_total_reads(self):
        with tempfile.TemporaryDirectory() as temporary:
            inputs, outputs = self.fixture(Path(temporary))
            path = outputs / "assignment_coverage.tsv"
            path.write_text(path.read_text().replace("\t6\n", "\t7\n"))
            with self.assertRaisesRegex(ValueError, "Coverage totals"):
                validate_taxonomy_tables(outputs, inputs, self.config)
