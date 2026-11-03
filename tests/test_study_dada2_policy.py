"""Independent DADA2 policy and input binding; synthetic mutations only."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/study_dada2"))
import execute


class StudyDada2PolicyTests(unittest.TestCase):
    def test_registered_policy_uses_terminal_inputs_and_original_parameters(self):
        path = ROOT / "config/prjna627078_dada2.yaml"
        policy, terminal, scope, manifest, processing = execute.load_policy(path)
        self.assertEqual(len(manifest), policy["expected_runs"])
        self.assertEqual(scope["processing_ceiling"], "qc")
        command = execute.dada_command(policy, terminal, scope)
        self.assertEqual(command[command.index("--input-dir") + 1], terminal["interim_dir"])
        self.assertEqual(command[command.index("--config") + 1], scope["project_config"])
        self.assertEqual(processing["dada2"]["filter"]["min_len"], 150)
        paths = execute.inputs(path, policy, terminal, scope, manifest)
        self.assertEqual(sum(p.name.endswith(".fastq.gz") for p in paths), 2 * len(manifest))
        self.assertIn(ROOT / "scripts/study_dada2/cohort.py", paths)
        self.assertIn(ROOT / "scripts/dada2/helpers.R", paths)

    def test_policy_rejects_cohort_changes_overrides_and_overlapping_outputs(self):
        policy = json.loads((ROOT / "config/prjna627078_dada2.yaml").read_text())
        mutations = [("expected_runs", 59), ("expected_batches", 3), ("threads", 4),
            ("sample_exclusions", True), ("parameter_overrides", True),
            ("require_convergence_before_limit", False), ("require_last_error_matrix_match", False),
            ("output_dir", "results/dada2/mota_gutierrez_2021_mexico/full"),
            ("audit_dir", policy["output_dir"] + "/audit"),
            ("intermediate_dir", "results/intermediate/dada2/other_study/terminal_3prime"),
            ("expected_batch_times", True)]
        for key, value in mutations:
            altered = copy.deepcopy(policy)
            altered[key] = value
            with self.subTest(key=key), tempfile.TemporaryDirectory() as folder:
                path = Path(folder) / "policy.json"
                path.write_text(json.dumps(altered))
                with self.assertRaises(ValueError):
                    execute.load_policy(path)
