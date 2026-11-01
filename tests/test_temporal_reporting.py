"""Synthetic temporal grouping, identity and missing-observation safeguards."""

import copy
import importlib.util
import json
import shutil
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
SPEC = importlib.util.spec_from_file_location("plot_temporal_diversity", ROOT / "scripts/reporting/plot_temporal_diversity.py")
REPORTING = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORTING)


class TemporalReportingTests(unittest.TestCase):
    def setUp(self):
        base = dict(study_id="study", run_accession="SRR1", fermentation_batch="study::batch", relative_time="0", fermentation_stage="early", sampling_stratum="1", observed_asvs="2", shannon_effective="1.1", library_reads="100")
        self.rows = [dict(base, run_accession="SRR2", sample_id="later", fermentation_hours="12"), dict(base, sample_id="start", fermentation_hours="0")]

    def test_sorting_preserves_values_and_sample_identity(self):
        result = REPORTING.prepare_plot_rows(self.rows, 2, 1)
        self.assertEqual([row["sample_id"] for row in result], ["start", "later"])
        self.assertEqual({row["sample_id"]: row for row in result}, {row["sample_id"]: row for row in self.rows})

    def test_duplicate_batch_hour_stratum_does_not_silently_pool(self):
        self.rows[0]["fermentation_hours"] = "0.0"
        with self.assertRaisesRegex(ValueError, "Duplicate batch/hour/stratum"):
            REPORTING.prepare_plot_rows(self.rows, 2, 1)

    def test_duplicate_sample_or_run_identity_fails(self):
        for field in ("sample_id", "run_accession"):
            altered = copy.deepcopy(self.rows)
            altered[0][field] = altered[1][field]
            with self.subTest(field=field), self.assertRaises(ValueError):
                REPORTING.prepare_plot_rows(altered, 2, 1)

    def test_strata_at_same_time_remain_separate(self):
        self.rows[0].update(fermentation_hours="0", sampling_stratum="2")
        result = REPORTING.prepare_plot_rows(self.rows, 2, 1)
        self.assertEqual(len(result), 2)
        self.assertEqual(REPORTING.table_summary(result, True)["series"], 2)

    def test_mixed_studies_missing_metadata_and_nonfinite_metrics_fail(self):
        for field, value in (("study_id", "other"), ("sampling_stratum", ""), ("fermentation_hours", "nan"), ("shannon_effective", "inf"), ("library_reads", "1.5")):
            altered = copy.deepcopy(self.rows)
            altered[0][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                REPORTING.prepare_plot_rows(altered, 2, 1)

    def test_profiles_reject_inference_aggregation_unsafe_paths_and_low_dpi(self):
        config = json.loads((ROOT / "config/temporal_reporting.yaml").read_text())
        REPORTING.temporal_profile(config, "full_study", ROOT)
        for field, value in (("inference", True), ("sample_aggregation", "batch"), ("break_at_missing_batch_times", False), ("x_field", "sample_index")):
            altered = copy.deepcopy(config)
            altered[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                REPORTING.temporal_profile(altered, "full_study", ROOT)
        config["figures"]["png_dpi"] = 72
        with self.assertRaises(ValueError):
            REPORTING.temporal_profile(config, "full_study", ROOT)
        config["figures"]["png_dpi"] = 300
        config["analyses"]["full_study"]["output_dir"] = "data/raw/wrong"
        with self.assertRaises(ValueError):
            REPORTING.temporal_profile(config, "full_study", ROOT)

    @unittest.skipUnless(shutil.which("Rscript"), "Requires the registered R environment")
    def test_r_segments_break_at_unobserved_batch_times(self):
        code = "source('scripts/reporting/plot_temporal_diversity.R'); p <- adjacent_pairs(c(0,12,36), c(0,12,24,36)); stopifnot(nrow(p)==1L, all(p[1,]==c(1,2)), nrow(adjacent_pairs(0, c(0,12)))==0L); failed <- inherits(try(adjacent_pairs(c(12,0), c(0,12)), silent=TRUE),'try-error'); stopifnot(failed)"
        result = subprocess.run(["Rscript", "-e", code], cwd=ROOT, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
