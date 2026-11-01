"""Cell comparison audits synthetic data; never rewrites classification tables."""

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/taxonomy"))
from audit_repetition import TABLES, compare_tables


class TaxonomyRepetitionTests(unittest.TestCase):
    def fixture(self, directory):
        for label in ("before", "after"):
            target = directory / label
            target.mkdir()
            for name in TABLES:
                (target / name).write_text("study_id\tasv_id\tGenus\nfixture\tasv1\tG1\nfixture\tasv2\tG2\n")
        return directory / "before", directory / "after"

    def test_cell_audit_separates_support_and_label_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            before, after = self.fixture(Path(temporary))
            target = after / "taxonomy_bootstraps.tsv"
            target.write_text(target.read_text().replace("G1", "G3"))
            audit = compare_tables(before, after)
            self.assertTrue(audit["taxonomy.tsv"]["identical"])
            self.assertEqual(audit["taxonomy_bootstraps.tsv"]["changed_cells"], 1)
            self.assertEqual(audit["taxonomy_bootstraps.tsv"]["changed_rows"], 1)

    def test_missing_row_fails_instead_of_becoming_zero_differences(self):
        with tempfile.TemporaryDirectory() as temporary:
            before, after = self.fixture(Path(temporary))
            (after / "taxonomy.tsv").write_text("study_id\tasv_id\tGenus\nfixture\tasv1\tG1\n")
            with self.assertRaisesRegex(ValueError, "row keys"):
                compare_tables(before, after)
