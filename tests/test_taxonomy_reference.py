"""Synthetic references test input validation; they are not scientific evidence."""

import copy
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.config import load_json_yaml
from cacao_inventory.taxonomy_reference import inspect_training_fasta, validate_taxonomy_config


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.config = load_json_yaml(ROOT / "config/taxonomy.yaml")

    def test_pinned_configuration_is_valid(self):
        validate_taxonomy_config(self.config, ROOT)

    def test_reference_path_cannot_escape_ignored_database_directory(self):
        for path in ("../outside.fa.gz", "references/reference.fa.gz"):
            config = copy.deepcopy(self.config)
            config["reference"]["path"] = path
            with self.assertRaises(ValueError):
                validate_taxonomy_config(config, ROOT)

    def test_invalid_threshold_and_automatic_exclusion_fail(self):
        for section, key, value in (("classification", "primary_min_boot", 101),
                                    ("screening", "exclude_flagged_asvs", True)):
            config = copy.deepcopy(self.config)
            config[section][key] = value
            with self.assertRaises(ValueError):
                validate_taxonomy_config(config, ROOT)

    def inspect(self, text):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fixture.fa.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(text)
            return inspect_training_fasta(path, 6)

    def test_multiline_fasta_and_incomplete_lineage_are_preserved(self):
        self.assertEqual(self.inspect(">Bacteria;P;C;O;F;G;\nACGT\nACGT\n>Bacteria;P;\nNN\n"),
                         {"reference_sequences": 2, "reference_bases": 10})

    def test_empty_sequence_invalid_bases_and_extra_rank_fail(self):
        for text in (">Bacteria;\n>Bacteria;\nACGT\n", ">Bacteria;\nACGZ\n", ">B;P;C;O;F;G;S;\nACGT\n"):
            with self.assertRaises(ValueError):
                self.inspect(text)
