"""Synthetic artifacts exercise cross-file checks; they are not biological evidence."""

from __future__ import annotations

import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/dada2"))
from validate_outputs import validate_run  # noqa: E402


def write_table(directory: Path, name: str, rows: list[dict[str, str]]) -> None:
    with (directory / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def checksum_outputs(directory: Path) -> None:
    rows = [{"path": path.name, "md5": hashlib.md5(path.read_bytes()).hexdigest()}
            for path in sorted(directory.iterdir()) if path.name not in {"SUCCESS", "output_checksums.tsv"}]
    write_table(directory, "output_checksums.tsv", rows)


class Dada2ArtifactTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)
        (self.directory / "SUCCESS").write_text("a" * 40 + "\n")
        write_table(self.directory, "sample_metadata.tsv", [{"study_id": "fixture", "sample_id": "s1", "run_accession": "SRR1"}])
        write_table(self.directory, "read_tracking.tsv", [{"study_id": "fixture", "sample_id": "s1", "run_accession": "SRR1",
                    "input": "10", "filtered": "9", "denoised_f": "9", "denoised_r": "8", "merged": "7", "nonchim": "6", "asv_count": "1"}])
        write_table(self.directory, "asv_counts.tsv", [{"study_id": "fixture", "sample_id": "s1", "fixture__ASV1": "6"}])
        write_table(self.directory, "asv_sequences.tsv", [{"study_id": "fixture", "asv_id": "fixture__ASV1", "sequence": "ACGT", "length": "4", "total_reads": "6"}])
        (self.directory / "asv_sequences.fasta").write_text(">fixture__ASV1\nACGT\n")
        write_table(self.directory, "summary.tsv", [{"study_id": "fixture", "samples": "1", "merged_asvs": "1", "nonchim_asvs": "1", "input_pairs": "10", "filtered_pairs": "9", "merged_pairs": "7", "nonchim_pairs": "6"}])
        write_table(self.directory, "sequence_length_distribution.tsv", [{"stage": stage, "length": "4", "asv_count": "1", "read_count": count} for stage, count in (("merged", "7"), ("nonchim", "6"))])
        write_table(self.directory, "error_learning.tsv", [{"read_direction": value, "converged_before_limit": "TRUE"} for value in ("F", "R")])
        (self.directory / "config_snapshot.yaml").write_text("fixture: true\n")
        write_table(self.directory, "software_versions.tsv", [{"software": "fixture", "version": "1"}])
        write_table(self.directory, "input_checksums.tsv", [{"path": "config_snapshot.yaml", "md5": hashlib.md5((self.directory / "config_snapshot.yaml").read_bytes()).hexdigest()}])
        checksum_outputs(self.directory)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_consistent_counts_and_sequences_pass(self) -> None:
        result = validate_run(self.directory)
        self.assertEqual(result["nonchim_pairs"], 6)
        self.assertTrue(result["error_models_converged_before_limit"])

    def test_changed_artifact_fails_checksum(self) -> None:
        with (self.directory / "asv_counts.tsv").open("a") as handle:
            handle.write("altered\n")
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            validate_run(self.directory)

    def test_rechecks_count_conservation_after_valid_checksums(self) -> None:
        write_table(self.directory, "asv_counts.tsv", [{"study_id": "fixture", "sample_id": "s1", "fixture__ASV1": "5"}])
        checksum_outputs(self.directory)
        with self.assertRaisesRegex(ValueError, "retained pairs"):
            validate_run(self.directory)

    def test_incomplete_run_cannot_pass(self) -> None:
        (self.directory / "SUCCESS").unlink()
        with self.assertRaises(FileNotFoundError):
            validate_run(self.directory)

    def test_length_histogram_must_match_sequences(self) -> None:
        write_table(self.directory, "sequence_length_distribution.tsv", [{"stage": "nonchim", "length": "5", "asv_count": "1", "read_count": "6"}])
        checksum_outputs(self.directory)
        with self.assertRaisesRegex(ValueError, "length distribution disagrees"):
            validate_run(self.directory)

    def test_original_inputs_can_be_verified(self) -> None:
        result = validate_run(self.directory, self.directory)
        self.assertEqual(result["inputs_checked"], 1)

    def test_stale_inputs_fail_even_with_updated_output_checksums(self) -> None:
        (self.directory / "config_snapshot.yaml").write_text("fixture: changed\n")
        checksum_outputs(self.directory)
        with self.assertRaisesRegex(ValueError, "Input checksum mismatch"):
            validate_run(self.directory, self.directory)

    def test_input_paths_cannot_escape_repository(self) -> None:
        write_table(self.directory, "input_checksums.tsv", [{"path": "../outside", "md5": "0" * 32}])
        checksum_outputs(self.directory)
        with self.assertRaisesRegex(ValueError, "input path escapes"):
            validate_run(self.directory, self.directory)


if __name__ == "__main__":
    unittest.main()
