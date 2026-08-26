from __future__ import annotations

import csv
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.qc import (  # noqa: E402
    build_raw_quality_rows,
    cutadapt_detection_counts,
    read_direction,
    validate_iupac_sequence,
)


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

    def test_cutadapt_primer_counts_are_cross_checked(self) -> None:
        report = {
            "read_counts": {"input": 10, "read1_with_adapter": 7},
            "adapters_read1": [{"total_matches": 7}],
        }
        self.assertEqual(cutadapt_detection_counts(report), (10, 7))
        report["adapters_read1"][0]["total_matches"] = 6
        with self.assertRaisesRegex(ValueError, "counts disagree"):
            cutadapt_detection_counts(report)

    def test_invalid_primer_symbols_fail(self) -> None:
        with self.assertRaisesRegex(ValueError, "invalid IUPAC"):
            validate_iupac_sequence("ACGTZ")


class PilotPrimerDetectionOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with (ROOT / "results" / "qc" / "pilot" / "primer_detection.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            cls.rows = list(csv.DictReader(handle, delimiter="\t"))

    def test_all_fastq_orientation_cases_are_recorded(self) -> None:
        self.assertEqual(len(self.rows), 96)
        read_files = {row["read_file"] for row in self.rows}
        self.assertEqual(len(read_files), 12)
        for read_file in read_files:
            self.assertEqual(
                sum(row["read_file"] == read_file for row in self.rows),
                8,
            )

    def test_detection_counts_are_bounded(self) -> None:
        for row in self.rows:
            examined = int(row["reads_examined"])
            matched = int(row["primer_matches"])
            self.assertGreater(examined, 0)
            self.assertGreaterEqual(matched, 0)
            self.assertLessEqual(matched, examined)

    def test_expected_full_constructs_are_absent(self) -> None:
        expected_constructs = [
            row
            for row in self.rows
            if row["expected_case"] == "true"
            and row["sequence_type"] == "full_construct"
        ]
        self.assertEqual(len(expected_constructs), 12)
        self.assertEqual(sum(int(row["primer_matches"]) for row in expected_constructs), 0)


if __name__ == "__main__":
    unittest.main()
