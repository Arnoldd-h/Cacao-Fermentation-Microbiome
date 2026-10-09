from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.design import build_analysis_units  # noqa: E402


class LongitudinalDesignTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = {"analysis_design": {
            "independent_unit": "fermentation_batch",
            "observation_fields": ["study_id", "fermentation_batch", "fermentation_hours"],
        }}
        self.studies = [{"study_id": "study", "bioproject": "PRJNA1", "include": "true",
                         "country": "country", "region_16s": "V4"}]
        self.samples = [{
            "study_id": "study", "bioproject": "PRJNA1", "sample_id": f"sample_{i}",
            "run_accession": f"SRR{i}", "analysis_include": "true",
            "fermentation_batch": "batch_1", "fermentation_hours": hour,
        } for i, hour in enumerate(("0", "0.0", "24", "24.0"), 1)]

    def test_subsamples_and_times_are_not_independent_batches(self) -> None:
        units, summaries = build_analysis_units(self.samples, self.studies, self.config)
        self.assertEqual(len(units), 4)
        self.assertEqual(summaries[0]["independent_batches"], "1")
        self.assertEqual(summaries[0]["batch_time_observations"], "2")
        self.assertEqual({row["subsamples_in_observation"] for row in units}, {"2"})

    def test_missing_batch_fails_instead_of_inventing_replicate(self) -> None:
        self.samples[0]["fermentation_batch"] = ""
        with self.assertRaisesRegex(ValueError, "Incomplete longitudinal unit"):
            build_analysis_units(self.samples, self.studies, self.config)

    def test_pending_samples_are_not_promoted(self) -> None:
        self.samples[0]["analysis_include"] = "pending"
        units, _ = build_analysis_units(self.samples, self.studies, self.config)
        self.assertEqual(len(units), 3)

    def test_duplicate_runs_fail(self) -> None:
        self.samples.append(self.samples[0].copy())
        with self.assertRaisesRegex(ValueError, "duplicate run"):
            build_analysis_units(self.samples, self.studies, self.config)

    def test_nonfinite_time_fails(self) -> None:
        self.samples[0]["fermentation_hours"] = "NaN"
        with self.assertRaisesRegex(ValueError, "Invalid fermentation_hours"):
            build_analysis_units(self.samples, self.studies, self.config)


if __name__ == "__main__":
    unittest.main()
