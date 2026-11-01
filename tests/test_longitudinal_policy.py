"""Registered policy checks; synthetic changes are never scientific evidence."""

import copy
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/longitudinal"))
from cacao_inventory.config import load_json_yaml
from longitudinal_policy import validate_longitudinal_policy


class LongitudinalPolicyTests(unittest.TestCase):
    def setUp(self):
        self.config = load_json_yaml(ROOT / "config/longitudinal.yaml")

    def test_registered_policy_and_six_scenarios(self):
        source, output = validate_longitudinal_policy(self.config, ROOT)
        self.assertEqual(len(self.config["scenarios"]), 6)
        self.assertTrue(output.is_relative_to(ROOT / "results/longitudinal"))
        self.assertTrue(source.is_file())

    def test_independence_and_pilot_inference_cannot_change_silently(self):
        for key, value in (("independent_unit", "sample_id"), ("pilot_inference", True), ("random_sampling", True)):
            config = copy.deepcopy(self.config)
            config["method"][key] = value
            with self.assertRaises(ValueError):
                validate_longitudinal_policy(config, ROOT)

    def test_invalid_numeric_settings_strata_and_paths_fail(self):
        for section, key, value in (("method", "alpha", True), ("method", "maximum_sign_patterns", 2),
                                    ("cohort", "expected_strata", ["1", "1"]),
                                    ("cohort", "minimum_timepoints_per_batch", 2)):
            config = copy.deepcopy(self.config)
            config[section][key] = value
            with self.assertRaises(ValueError):
                validate_longitudinal_policy(config, ROOT)
        config = copy.deepcopy(self.config)
        config["output_dir"] = "results/taxonomy/pilot"
        with self.assertRaises(ValueError):
            validate_longitudinal_policy(config, ROOT)

    def test_unknown_stratum_and_nonpositive_pseudocount_fail(self):
        for key, value in (("stratum", "unknown"), ("pseudocount", 0)):
            config = copy.deepcopy(self.config)
            config["scenarios"][-1][key] = value
            with self.assertRaises(ValueError):
                validate_longitudinal_policy(config, ROOT)
