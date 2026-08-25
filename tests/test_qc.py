from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.qc import build_raw_quality_rows, read_direction  # noqa: E402


def validation_rows() -> list[dict[str, str]]:
    forward = {
        "study_id": "study",
        "bioproject": "PRJNA1",
        "sample_id": "sample",
        "run_accession": "SRR1",
        "read_file": "SRR1_1.fastq.gz",
        "read_count": "10",
        "base_count": "1000",
        "minimum_read_length": "100",
        "maximum_read_length": "100",
        "status": "valid",
    }
    reverse = dict(forward)
    reverse["read_file"] = "SRR1_2.fastq.gz"
    return [forward, reverse]


def multiqc_rows() -> list[dict[str, str]]:
    forward = {
        "Sample": "SRR1_1",
        "Filename": "SRR1_1.fastq.gz",
        "Total Sequences": "10.0",
        "Sequences flagged as poor quality": "0.0",
        "Sequence length": "100",
        "%GC": "50.0",
        "avg_sequence_length": "100.0",
        "median_sequence_length": "100",
        "basic_statistics": "pass",
        "per_base_sequence_quality": "pass",
        "per_sequence_quality_scores": "pass",
        "per_base_sequence_content": "warn",
        "per_sequence_gc_content": "pass",
        "per_base_n_content": "pass",
        "sequence_length_distribution": "pass",
        "sequence_duplication_levels": "warn",
        "overrepresented_sequences": "warn",
        "adapter_content": "pass",
    }
    reverse = dict(forward)
    reverse["Sample"] = "SRR1_2"
    reverse["Filename"] = "SRR1_2.fastq.gz"
    return [forward, reverse]


class RawQualityTests(unittest.TestCase):
    def test_fastqc_rows_join_validated_fastq(self) -> None:
        observed = build_raw_quality_rows(multiqc_rows(), validation_rows())
        self.assertEqual(observed[0]["read_direction"], "R1")
        self.assertEqual(observed[0]["total_sequences"], "10")
        self.assertEqual(observed[0]["total_bases"], "1000")
        self.assertEqual(observed[1]["read_direction"], "R2")

    def test_fastqc_count_mismatch_fails(self) -> None:
        fastqc = multiqc_rows()
        fastqc[0]["Total Sequences"] = "9.0"
        with self.assertRaisesRegex(ValueError, "read count mismatch"):
            build_raw_quality_rows(fastqc, validation_rows())

    def test_incomplete_pair_fails(self) -> None:
        with self.assertRaisesRegex(ValueError, "directions are incomplete"):
            build_raw_quality_rows(multiqc_rows()[:1], validation_rows()[:1])

    def test_only_paired_ena_filenames_are_accepted(self) -> None:
        self.assertEqual(read_direction("run_2.fastq.gz"), "R2")
        with self.assertRaisesRegex(ValueError, "paired read direction"):
            read_direction("single.fastq.gz")


if __name__ == "__main__":
    unittest.main()
