from __future__ import annotations

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

from cacao_inventory.pilot_workflow import fastq_path, load_pilot_manifest, manifest_row, processing_configuration


class PilotWorkflowTests(unittest.TestCase):
    def test_processing_uses_manifest_bioproject(self) -> None:
        config = {"amplicon_processing": {
            "PRJNA1": {"study_id": "one", "primers": "first"},
            "PRJNA2": {"study_id": "two", "primers": "second"},
        }}
        observed = processing_configuration(config, {"study_id": "two", "bioproject": "PRJNA2"})
        self.assertEqual(observed["primers"], "second")
        with self.assertRaisesRegex(ValueError, "study_id disagrees"):
            processing_configuration(config, {"study_id": "one", "bioproject": "PRJNA2"})

    def test_manifest_paths_are_study_scoped(self) -> None:
        row = {"study_id": "study", "run_accession": "SRR1"}
        self.assertEqual(fastq_path(row, "2", "trimmed"), "data/interim/study/pilot/SRR1/SRR1_2.fastq.gz")
        with self.assertRaisesRegex(ValueError, "Study wildcard disagrees"):
            manifest_row([row], "SRR1", "other")

    def test_real_manifest_is_valid_and_paired(self) -> None:
        rows = load_pilot_manifest(ROOT / "metadata/pilot_manifest.tsv")
        self.assertEqual(len(rows), 6)


@unittest.skipUnless(shutil.which("snakemake"), "Snakemake integration checks require the declared Linux environment")
class WorkflowBootstrapTests(unittest.TestCase):
    def copy_sources(self, destination: Path) -> None:
        for directory in ("workflow", "scripts", "python", "config", "environment"):
            shutil.copytree(ROOT / directory, destination / directory, ignore=shutil.ignore_patterns("__pycache__"))
        (destination / "metadata").mkdir()
        for path in (ROOT / "metadata").glob("*.tsv"):
            shutil.copy2(path, destination / "metadata" / path.name)
        (destination / "results/tables").mkdir(parents=True)
        shutil.copy2(ROOT / "results/tables/pilot_dataset_selection.tsv", destination / "results/tables/pilot_dataset_selection.tsv")
        (destination / "metadata/.inventory.ok").write_text("validated test copy\n", encoding="utf-8")
        # Copied metadata are intentional inputs; no network/source refresh is
        # required to exercise the sequencing DAG in this isolated copy.
        newest = time.time() + 10
        for directory in (destination / "metadata", destination / "results/tables"):
            for path in directory.iterdir():
                os.utime(path, (newest, newest))

    def dry_run(self, root: Path, target: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["snakemake", "--snakefile", "workflow/Snakefile", "--cores", "1", "--dry-run", target],
            cwd=root, capture_output=True, text=True, timeout=90,
        )

    def test_absent_manifest_reaches_checkpoint_without_parse_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.copy_sources(root)
            (root / "metadata/pilot_manifest.tsv").unlink()
            result = self.dry_run(root, "pilot_post_trim_qc")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("build_pilot_manifest", result.stdout)

    def test_absent_raw_fastq_have_download_producers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.copy_sources(root)
            result = self.dry_run(root, "pilot_post_trim_qc")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("download_pilot_read", result.stdout)

    def test_existing_raw_is_not_replaced_when_metadata_is_newer(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.copy_sources(root)
            row = load_pilot_manifest(root / "metadata/pilot_manifest.tsv")[0]
            target = fastq_path(row, "1")
            raw = root / target
            raw.parent.mkdir(parents=True)
            raw.write_bytes(b"synthetic fixture for DAG planning only")
            older = time.time() - 100
            os.utime(raw, (older, older))
            raw.chmod(0o444)
            result = self.dry_run(root, target)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("Nothing to be done", result.stdout)
            self.assertEqual(raw.read_bytes(), b"synthetic fixture for DAG planning only")


if __name__ == "__main__":
    unittest.main()
