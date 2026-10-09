"""Full-study reconstruction and immutable raw-file regression checks."""

import base64
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
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import read_inventory, write_tsv_atomic
from cacao_inventory.study_scope import MANIFEST_COLUMNS, build_study_manifest


@unittest.skipUnless(shutil.which("snakemake"), "Requires the declared Linux environment")
class StudyWorkflowTests(unittest.TestCase):
    def copy_sources(self, root):
        for directory in ("workflow", "scripts", "python", "config", "environment"):
            shutil.copytree(ROOT / directory, root / directory, ignore=shutil.ignore_patterns("__pycache__"))
        (root / "metadata").mkdir()
        for path in (ROOT / "metadata").glob("*.tsv"):
            shutil.copy2(path, root / "metadata" / path.name)

    def dry_run(self, root, target="all"):
        return subprocess.run(["snakemake", "--snakefile", "workflow/study.Snakefile", "--cores", "1", "--dry-run", target], cwd=root, capture_output=True, text=True, timeout=90)

    def test_absent_manifest_reaches_source_checkpoint_and_reference(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.copy_sources(root)
            result = self.dry_run(root)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("study_manifest", result.stdout)
            self.assertIn("download_taxonomy_database", result.stdout)

    def test_shared_reference_keeps_existing_pilot_scheduler_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.copy_sources(root)
            reference = load_json_yaml(root / "config/taxonomy.yaml")["reference"]
            for name in ("path", "provenance_path"):
                relative = reference[name]
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"synthetic scheduler fixture only")
                encoded = base64.urlsafe_b64encode(relative.encode()).decode()
                history = ROOT / ".snakemake/metadata" / encoded
                if not history.exists():
                    self.skipTest("Requires completed pilot reference scheduler history")
                destination = root / ".snakemake/metadata" / encoded
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(history, destination)
            (root / reference["path"]).chmod(0o444)
            result = self.dry_run(root, reference["provenance_path"])
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("Nothing to be done", result.stdout)

    def test_registered_config_refresh_does_not_replace_existing_raw(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.copy_sources(root)
            cfg = load_json_yaml(root / "config/full_study.yaml")
            inventory = read_inventory(root / "metadata")
            manifest = build_study_manifest(cfg, inventory["runs"], inventory["samples"])
            write_tsv_atomic(root / cfg["manifest"], manifest, MANIFEST_COLUMNS)
            for name in ("source_audit.json", "resource_plan.json", "fastq_stream_log.tsv", "provenance.json", "input_checksums.json", "config_snapshot.yaml", "SUCCESS"):
                path = root / cfg["plan_dir"] / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("synthetic DAG fixture only\n")
            row = manifest[0]
            target = "data/raw/" + cfg["study_id"] + "/" + row["run_accession"] + "/" + row["run_accession"] + "_1.fastq.gz"
            raw = root / target
            raw.parent.mkdir(parents=True, exist_ok=True)
            raw.write_bytes(b"synthetic DAG fixture only")
            older = time.time() - 100
            os.utime(raw, (older, older))
            raw.chmod(0o444)
            result = self.dry_run(root, target)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertNotIn("rule download_read:", result.stdout)


if __name__ == "__main__":
    unittest.main()
