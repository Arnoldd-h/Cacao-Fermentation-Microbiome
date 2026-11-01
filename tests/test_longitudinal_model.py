"""Analytical/synthetic longitudinal checks; no fixture is biological evidence."""

import copy
import json
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/longitudinal"))
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from longitudinal_model import (IDENTITY, adjust_by, aggregate_observations, build_tables, clr_vector,
                                exact_sign_test, fit_multivariate_ols, mean_vector, table_columns, validate_cohort)
from validate_longitudinal import compare_table


def synthetic_cohort(times=(3, 3, 3, 3), multipliers=(1, 1, 1, 1)):
    config = copy.deepcopy(load_json_yaml(ROOT / "config/longitudinal.yaml"))
    config["cohort"].update(study_id="synthetic", expected_samples=2 * sum(times), expected_batches=len(times), expected_batch_times=sum(times))
    counts, metadata = [], []
    for batch_index, n_times in enumerate(times):
        batch = f"synthetic::batch_{batch_index}"
        for time in range(n_times):
            for stratum in ("1", "2"):
                sample = f"batch_{batch_index}_time_{time}_stratum_{stratum}"
                counts.append({"study_id": "synthetic", "sample_id": sample,
                    "ASV1": str(2 ** (multipliers[batch_index] * time + 1) - 1), "ASV2": "1"})
                metadata.append(dict(zip(IDENTITY, ["synthetic", sample, f"synthetic_run_{sample}", batch,
                    str(12 * time), str(time / (n_times - 1)), str(12 * (n_times - 1)), stratum])))
    return config, counts, ["ASV1", "ASV2"], metadata


class LongitudinalModelTests(unittest.TestCase):
    def test_ols_has_intercept_and_exact_linear_solution(self):
        slope, intercept, rss, total, r2 = fit_multivariate_ols([0, 0.5, 1], [[10, -10], [11, -11], [12, -12]])
        self.assertEqual(slope, [2, -2])
        self.assertEqual(intercept, [10, -10])
        self.assertEqual(rss, 0)
        self.assertEqual(total, 4)
        self.assertEqual(r2, 1)

    def test_clr_center_equal_strata_does_not_pool_library_depths(self):
        observed = mean_vector([clr_vector([90, 10], 1), clr_vector([1, 9], 1)])
        expected = math.log((91 * 2) / (11 * 10)) / 4
        self.assertAlmostEqual(observed[0], expected)
        self.assertAlmostEqual(sum(observed), 0)
        self.assertNotAlmostEqual(observed[0], clr_vector([91, 19], 1)[0])

    def test_analytic_clr_slope_and_global_exact_resolution(self):
        config, counts, features, metadata = synthetic_cohort()
        tables = build_tables(counts, features, metadata, config)
        primary = tables["global_tests.tsv"][0]
        self.assertAlmostEqual(primary["norm_mean_slope"], math.sqrt(2) * math.log(2))
        self.assertEqual(primary["exact_p_value"], 0.125)
        self.assertEqual(primary["sign_patterns"], 16)
        self.assertEqual(primary["n_samples"], 24)
        for row in tables["slopes.tsv"][:4]:
            self.assertAlmostEqual(row["ASV1"], math.log(2))
            self.assertAlmostEqual(row["ASV2"], -math.log(2))

    def test_equal_batch_weights_with_unequal_timepoints(self):
        config, counts, features, metadata = synthetic_cohort((3, 4, 5, 6), (1, 2, 3, 4))
        primary = build_tables(counts, features, metadata, config)["global_tests.tsv"][0]
        # Analytic slopes are (1,3,6,10)*log(2), regardless of time counts.
        expected = 5 * math.sqrt(2) * math.log(2)
        self.assertAlmostEqual(primary["norm_mean_slope"], expected)
        weighted_wrong = (3 * 1 + 4 * 3 + 5 * 6 + 6 * 10) / 18 * math.sqrt(2) * math.log(2)
        self.assertNotAlmostEqual(primary["norm_mean_slope"], weighted_wrong)

    def test_incomplete_observation_audited_and_retained_in_sensitivity(self):
        config, counts, features, metadata = synthetic_cohort((4, 4, 4, 4))
        missing = "batch_0_time_1_stratum_2"
        counts = [row for row in counts if row["sample_id"] != missing]
        metadata = [row for row in metadata if row["sample_id"] != missing]
        config["cohort"]["expected_samples"] -= 1
        tables = build_tables(counts, features, metadata, config)
        primary, available = tables["global_tests.tsv"][0], tables["global_tests.tsv"][3]
        self.assertEqual((primary["n_observations"], primary["n_samples"]), (15, 30))
        self.assertEqual((available["n_observations"], available["n_samples"]), (16, 31))
        rejected = [row for row in tables["sample_membership.tsv"] if row["scenario_id"] == config["primary_scenario"] and row["included"] == "FALSE"]
        self.assertEqual(len(rejected), 1)
        self.assertEqual(rejected[0]["reason"], "incomplete_strata_observation")

    def test_invalid_duration_duplicate_stratum_and_empty_library_fail(self):
        for mutation in ("duration", "duplicate", "empty"):
            config, counts, features, metadata = synthetic_cohort()
            if mutation == "duration":
                metadata[0]["fermentation_duration_hours"] = "0"
            elif mutation == "duplicate":
                metadata[1]["sampling_stratum"] = "1"
            else:
                counts[0].update(ASV1="0", ASV2="0")
            with self.assertRaises(ValueError):
                validate_cohort(counts, features, metadata, config)

    def test_missing_baseline_or_too_few_times_fails_without_relaxation(self):
        config, counts, features, metadata = synthetic_cohort()
        metadata = [row for row in metadata if row["sample_id"] != "batch_0_time_0_stratum_2"]
        counts = [row for row in counts if row["sample_id"] != "batch_0_time_0_stratum_2"]
        config["cohort"]["expected_samples"] -= 1
        with self.assertRaisesRegex(ValueError, "baseline and enough"):
            aggregate_observations(counts, features, metadata, config, config["scenarios"][0])

    def test_sign_test_zero_or_opposite_slopes_and_by(self):
        _, p, patterns = exact_sign_test([[0, 0]] * 4, 1e-12, 16)
        self.assertEqual((p, len(patterns)), (1, 16))
        self.assertEqual(exact_sign_test([[1, -1], [-1, 1], [1, -1], [-1, 1]], 1e-12, 16)[1], 1)
        expected = [0.055, 0.07333333333333333, 0.07333333333333333]
        for observed, target in zip(adjust_by([0.01, 0.04, 0.03]), expected):
            self.assertAlmostEqual(observed, target)

    def test_table_validator_catches_changed_coefficient(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "slopes.tsv"
            write_tsv_atomic(path, [{"id": "synthetic", "slope": 4}], ["id", "slope"])
            with self.assertRaisesRegex(ValueError, "formula mismatch"):
                compare_table(path, [{"id": "synthetic", "slope": 2}], ["id", "slope"], 1e-8)

    @unittest.skipUnless(shutil.which("Rscript"), "Requires the declared R environment")
    def test_independent_r_qr_signs_and_by_then_rejects_tampered_slope(self):
        config, counts, features, metadata = synthetic_cohort((3, 4, 5, 6), (1, 2, 3, 4))
        tables = build_tables(counts, features, metadata, config)
        with tempfile.TemporaryDirectory(prefix="synthetic_longitudinal_") as temporary:
            root = Path(temporary)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config) + "\n")
            write_tsv_atomic(root / "counts.tsv", counts, ["study_id", "sample_id", *features])
            write_tsv_atomic(root / "metadata.tsv", metadata, IDENTITY)
            for name, columns in table_columns(features).items():
                write_tsv_atomic(root / name, tables[name], columns)
            command = ["Rscript", "scripts/longitudinal/validate_model.R", "--config", str(config_path), "--counts", str(root / "counts.tsv"),
                       "--metadata", str(root / "metadata.tsv"), "--output-dir", str(root)]
            completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=90)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            tables["slopes.tsv"][0]["ASV1"] += 0.1
            write_tsv_atomic(root / "slopes.tsv", tables["slopes.tsv"], table_columns(features)["slopes.tsv"])
            failed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=90)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("QR slopes", failed.stderr)
