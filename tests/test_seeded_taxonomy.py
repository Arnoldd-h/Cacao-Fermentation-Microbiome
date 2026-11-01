"""Repeat native classification on synthetic ties, never on scientific fixtures."""

import copy
import hashlib
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.config import load_json_yaml
sys.path.insert(0, str(ROOT / "scripts/taxonomy"))
from seeded_policy import validate_native_rng


class SeededTaxonomyTests(unittest.TestCase):
    def test_upstream_source_is_exactly_pinned(self):
        config = load_json_yaml(ROOT / "config/taxonomy.yaml")
        expected = config["classification"]["native_rng"]["upstream_sha256"]
        self.assertEqual(hashlib.sha256((ROOT / "scripts/taxonomy/native/taxonomy_upstream.cpp").read_bytes()).hexdigest(), expected)

    def test_unregistered_native_adapter_fails(self):
        config = copy.deepcopy(load_json_yaml(ROOT / "config/taxonomy.yaml"))
        config["classification"]["native_rng"]["adapter"] = "unregistered"
        with self.assertRaisesRegex(ValueError, "Unregistered"):
            validate_native_rng(config)

    @unittest.skipUnless(shutil.which("Rscript"), "Requires the declared R/Rcpp compiler environment")
    def test_exact_repetition_in_separate_processes_and_thread_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            results = []
            for iteration, threads in enumerate((1, 2, 2)):
                output = Path(directory) / f"attempt_{iteration}.rds"
                completed = subprocess.run(
                    ["Rscript", "tests/taxonomy/test_seeded_classifier.R", str(output), str(threads)],
                    cwd=ROOT, text=True, capture_output=True, timeout=240)
                self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
                results.append(output.read_bytes())
            self.assertEqual(results[0], results[1])
            self.assertEqual(results[1], results[2])
