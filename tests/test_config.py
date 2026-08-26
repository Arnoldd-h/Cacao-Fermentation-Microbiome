from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import (  # noqa: E402
    ConfigurationError,
    assign_temporal_stage,
    load_project_configuration,
    validate_temporal_stages,
)
from cacao_inventory.qc import reverse_complement_iupac, validate_iupac_sequence  # noqa: E402


class TemporalStageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config, cls.candidates = load_project_configuration(
            ROOT / "config" / "config.yaml", ROOT / "config" / "datasets.yaml"
        )
        cls.stages = cls.config["temporal_stages"]

    def test_boundaries_are_right_closed(self) -> None:
        self.assertEqual(assign_temporal_stage(0.0, self.stages), "early")
        self.assertEqual(assign_temporal_stage(0.33, self.stages), "early")
        self.assertEqual(assign_temporal_stage(0.330001, self.stages), "mid")
        self.assertEqual(assign_temporal_stage(0.66, self.stages), "mid")
        self.assertEqual(assign_temporal_stage(0.660001, self.stages), "late")
        self.assertEqual(assign_temporal_stage(1.0, self.stages), "late")

    def test_missing_relative_time_has_no_stage(self) -> None:
        self.assertEqual(assign_temporal_stage(None, self.stages), "")

    def test_out_of_range_time_fails(self) -> None:
        with self.assertRaises(ValueError):
            assign_temporal_stage(1.01, self.stages)

    def test_non_contiguous_stages_fail(self) -> None:
        with self.assertRaises(ConfigurationError):
            validate_temporal_stages(
                {"early": [0, 0.2], "mid": [0.3, 0.6], "late": [0.6, 1]}
            )

    def test_candidate_accessions_are_unique(self) -> None:
        accessions = [candidate["bioproject"] for candidate in self.candidates]
        self.assertEqual(len(accessions), len(set(accessions)))

    def test_pilot_primer_configuration_is_exact_and_valid(self) -> None:
        pilot = self.config["amplicon_processing"]["PRJNA492720"]
        self.assertEqual(pilot["marker"], "16S rRNA")
        self.assertEqual(pilot["region"], "V4")
        expected = {
            "forward": "GTGCCAGCMGCCGCGGTAA",
            "reverse": "GGACTACHVGGGTWTCTAAT",
        }
        for role, sequence in expected.items():
            primer = pilot["primers"][role]
            self.assertEqual(primer["gene_specific_sequence"], sequence)
            self.assertEqual(validate_iupac_sequence(sequence), sequence)
            self.assertTrue(primer["full_construct_sequence"].endswith(sequence))

    def test_iupac_reverse_complements_are_stable(self) -> None:
        sequence = "GGACTACHVGGGTWTCTAAT"
        reverse = reverse_complement_iupac(sequence)
        self.assertEqual(reverse, "ATTAGAWACCCBDGTAGTCC")
        self.assertEqual(reverse_complement_iupac(reverse), sequence)


if __name__ == "__main__":
    unittest.main()
