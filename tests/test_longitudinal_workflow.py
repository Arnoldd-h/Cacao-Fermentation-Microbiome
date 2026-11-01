"""Synthetic scheduler fixtures only; touch never executes a scientific model."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/longitudinal"))
from cacao_inventory.config import load_json_yaml
from validate_longitudinal import PRODUCTS


@unittest.skipUnless(shutil.which("snakemake"), "Requires the declared workflow environment")
class LongitudinalWorkflowTests(unittest.TestCase):
    def fixture(self, root):
        for directory in ("workflow", "scripts", "python", "config", "environment"):
            shutil.copytree(ROOT / directory, root / directory, ignore=shutil.ignore_patterns("__pycache__"))
        config = load_json_yaml(root / "config/longitudinal.yaml")
        source = load_json_yaml(root / config["source_diversity_config"])["paths"]
        for directory, names in ((source["bacterial_dir"], ("asv_counts.tsv", "sample_metadata.tsv", "taxonomy.tsv")),
                                  (source["filtering_dir"], ("SUCCESS", "validation.json", "provenance.json"))):
            target = root / directory
            target.mkdir(parents=True, exist_ok=True)
            for name in names:
                (target / name).write_text("synthetic scheduler fixture only\n")
        return root / "config/longitudinal.yaml"

    def run_workflow(self, root, *args):
        return subprocess.run(["snakemake", "--snakefile", "workflow/longitudinal.Snakefile", "--cores", "1", *args],
                              cwd=root, capture_output=True, text=True, timeout=90)

    def test_only_downstream_model_jobs_are_scheduled(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.fixture(root)
            completed = self.run_workflow(root, "--dry-run")
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertIn("run_longitudinal", completed.stdout)
            for unexpected in ("run_taxonomy", "run_dada2", "download_read", "fastqc"):
                self.assertNotIn(f"rule {unexpected}:", completed.stdout)

    def test_config_content_change_invalidates_even_with_older_timestamp(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = self.fixture(root)
            config = load_json_yaml(path)
            output = root / config["output_dir"]
            output.mkdir(parents=True, exist_ok=True)
            for name in [*PRODUCTS, "config_snapshot.yaml", "input_checksums.json", "provenance.json", "SUCCESS", "validation.json"]:
                (output / name).write_text("synthetic scheduler fixture only\n")
            touched = self.run_workflow(root, "--touch", "--forceall")
            self.assertEqual(touched.returncode, 0, touched.stdout + touched.stderr)
            stable = self.run_workflow(root, "--dry-run")
            self.assertIn("Nothing to be done", stable.stdout)
            config = load_json_yaml(path)
            config["scenarios"][0]["pseudocount"] = 0.75
            path.write_text(json.dumps(config) + "\n")
            older = time.time() - 100
            os.utime(path, (older, older))
            changed = self.run_workflow(root, "--dry-run")
            self.assertEqual(changed.returncode, 0, changed.stdout + changed.stderr)
            self.assertIn("run_longitudinal", changed.stdout)
            self.assertIn("Params have changed", changed.stdout)
