from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.pilot import build_pilot_manifest  # noqa: E402
from cacao_inventory.schema import RUN_COLUMNS, SAMPLE_COLUMNS  # noqa: E402


RULES = {
    "samples_per_stage": 2,
    "required_stages": ["early", "mid", "late"],
    "required_layout": "PAIRED",
    "required_fastq_files": 2,
}


def blank(columns: list[str]) -> dict[str, str]:
    return {column: "" for column in columns}


def fixtures() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    runs: list[dict[str, str]] = []
    samples: list[dict[str, str]] = []
    for stage_index, stage in enumerate(RULES["required_stages"]):
        for batch_index in range(2):
            accession = f"SRR{stage_index}{batch_index}000001"
            run = blank(RUN_COLUMNS)
            run.update(
                {
                    "study_id": "pilot_study",
                    "bioproject": "PRJNA123456",
                    "run_accession": accession,
                    "library_layout": "PAIRED",
                    "sequencing_platform": "ILLUMINA",
                    "fastq_ftp": f"example.org/{accession}_1.fastq.gz;example.org/{accession}_2.fastq.gz",
                    "fastq_md5": "a" * 32 + ";" + "b" * 32,
                    "fastq_bytes": "100;200",
                    "analysis_include": "true",
                }
            )
            sample = blank(SAMPLE_COLUMNS)
            sample.update(
                {
                    "study_id": "pilot_study",
                    "bioproject": "PRJNA123456",
                    "sample_id": f"sample_{stage}_{batch_index}",
                    "run_accession": accession,
                    "fermentation_batch": f"pilot_study::batch_{batch_index}",
                    "fermentation_hours": str(stage_index * 48),
                    "relative_time": str(stage_index / 2),
                    "fermentation_stage": stage,
                    "sampling_stratum": str(batch_index + 1),
                    "analysis_include": "true",
                }
            )
            runs.append(run)
            samples.append(sample)
    return runs, samples


class PilotManifestTests(unittest.TestCase):
    def test_balanced_manifest_has_two_independent_batches_per_stage(self) -> None:
        runs, samples = fixtures()
        rows = build_pilot_manifest("pilot_study", runs, samples, RULES)

        self.assertEqual(len(rows), 6)
        self.assertEqual(sum(int(row["estimated_bytes_total"]) for row in rows), 1800)
        for stage in RULES["required_stages"]:
            stage_rows = [row for row in rows if row["fermentation_stage"] == stage]
            self.assertEqual(len(stage_rows), 2)
            self.assertEqual(len({row["fermentation_batch"] for row in stage_rows}), 2)
        self.assertTrue(rows[0]["fastq_ftp"].startswith("https://"))

    def test_incomplete_fastq_metadata_is_rejected(self) -> None:
        runs, samples = fixtures()
        runs[0]["fastq_md5"] = "a" * 32

        with self.assertRaisesRegex(ValueError, "Stage early has 1 eligible batches"):
            build_pilot_manifest("pilot_study", runs, samples, RULES)

    def test_malformed_md5_is_rejected(self) -> None:
        runs, samples = fixtures()
        runs[0]["fastq_md5"] = "not-an-md5;" + "b" * 32

        with self.assertRaisesRegex(ValueError, "Stage early has 1 eligible batches"):
            build_pilot_manifest("pilot_study", runs, samples, RULES)


if __name__ == "__main__":
    unittest.main()
