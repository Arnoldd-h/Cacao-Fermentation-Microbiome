from __future__ import annotations

import os
import gzip
import hashlib
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
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.schema import PILOT_MANIFEST_COLUMNS


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

    def test_manifest_timestamp_does_not_invalidate_but_changed_content_does(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.copy_sources(root)
            # A separate temporary repository and explicitly synthetic FASTQ
            # exercise real validation with no network and no project data edits.
            subprocess.run(["git", "init", "--quiet", str(root)], check=True)
            subprocess.run(["git", "-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "--allow-empty", "--quiet", "-m", "synthetic workflow fixture"], cwd=root, check=True)
            manifest = root / "metadata/pilot_manifest.tsv"
            rows = load_pilot_manifest(manifest)
            content = gzip.compress(b"@synthetic_fixture\nACGT\n+\nIIII\n", mtime=0)
            checksum = hashlib.md5(content, usedforsecurity=False).hexdigest()
            for row in rows:
                for direction in ("1", "2"):
                    path = root / fastq_path(row, direction)
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(content)
                row["fastq_ftp"] = ";".join(f"https://example.invalid/{row['run_accession']}_{direction}.fastq.gz" for direction in ("1", "2"))
                row["fastq_bytes"] = f"{len(content)};{len(content)}"
                row["fastq_md5"] = f"{checksum};{checksum}"
                row["estimated_bytes_total"] = str(2 * len(content))
            write_tsv_atomic(manifest, rows, PILOT_MANIFEST_COLUMNS)
            old = time.time() - 200
            for directory in ("workflow", "scripts", "python", "config", "environment"):
                for path in (root / directory).rglob("*"):
                    if path.is_file():
                        os.utime(path, (old, old))
            for directory in ("metadata", "results/tables", "data"):
                for path in (root / directory).rglob("*"):
                    if path.is_file():
                        os.utime(path, (old + 100, old + 100))
            result = subprocess.run(
                ["snakemake", "--snakefile", "workflow/Snakefile", "--cores", "1", "download_pilot_fastq"],
                cwd=root, capture_output=True, text=True, timeout=90,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            report = root / "results/qc/pilot_download_validation.tsv"
            self.assertIn("reused_valid", report.read_text(encoding="utf-8"))
            new_time = report.stat().st_mtime + 1
            os.utime(manifest, (new_time, new_time))
            unchanged = self.dry_run(root, "download_pilot_fastq")
            self.assertEqual(unchanged.returncode, 0, unchanged.stdout + unchanged.stderr)
            self.assertIn("Nothing to be done", unchanged.stdout)
            rows[0]["sample_id"] += "_changed_fixture"
            write_tsv_atomic(manifest, rows, PILOT_MANIFEST_COLUMNS)
            changed = self.dry_run(root, "download_pilot_fastq")
            self.assertEqual(changed.returncode, 0, changed.stdout + changed.stderr)
            self.assertIn("Params have changed", changed.stdout)
            self.assertNotIn("rule download_pilot_read:", changed.stdout)


if __name__ == "__main__":
    unittest.main()
