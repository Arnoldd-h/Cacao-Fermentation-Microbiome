"""Cell comparison audits synthetic data; never rewrites classification tables."""

import json
import subprocess
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

    def test_equal_tables_with_different_inputs_cannot_pass_exact_repetition(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            before, after = self.fixture(root)
            for index, directory in enumerate((before, after)):
                (directory / "SUCCESS").write_text("a" * 40 + "\n")
                (directory / "config_snapshot.yaml").write_text("{}\n")
                (directory / "provenance.json").write_text(json.dumps({"effective_threads": index + 1}))
                (directory / "input_checksums.json").write_text(json.dumps([
                    {"path": "synthetic_input.tsv", "sha256": str(index) * 64}]))
            completed = subprocess.run([sys.executable, str(ROOT / "scripts/taxonomy/audit_repetition.py"),
                "--baseline", str(before), "--repeat", str(after), "--output", str(root / "audit.json"),
                "--require-identical"], capture_output=True, text=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("inputs/configuration differ", completed.stderr)
