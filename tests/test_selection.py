from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.schema import RUN_COLUMNS, SAMPLE_COLUMNS, STUDY_COLUMNS  # noqa: E402
from cacao_inventory.selection import build_pilot_selection  # noqa: E402


RULES = {
    "minimum_timepoints": 3,
    "required_temporal_stages": ["early", "mid", "late"],
    "metadata_quality_order": ["high", "medium", "low"],
    "ranking_priority": [
        "metadata_quality",
        "primers_known",
        "paired_end",
        "unique_timepoints",
        "fermentation_batches",
        "selected_fastq_bytes",
    ],
}


def blank(columns: list[str]) -> dict[str, str]:
    return {column: "" for column in columns}


def fixture_study(study_id: str, bioproject: str) -> dict[str, str]:
    row = blank(STUDY_COLUMNS)
    row.update(
        {
            "study_id": study_id,
            "bioproject": bioproject,
            "doi": "10.0000/example",
            "marker": "16S rRNA",
            "forward_primer": "515F",
            "reverse_primer": "806R",
            "paired_end": "true",
            "raw_data_available": "true",
            "metadata_quality": "high",
            "include": "true",
        }
    )
    return row


def fixture_rows(study_id: str, bioproject: str, hours: list[int]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    stages = ["early", "mid", "late"]
    runs: list[dict[str, str]] = []
    samples: list[dict[str, str]] = []
    for index, hour in enumerate(hours):
        accession = f"SRR{study_id[-1]}{index:06d}"
        run = blank(RUN_COLUMNS)
        run.update(
            {
                "study_id": study_id,
                "bioproject": bioproject,
                "run_accession": accession,
                "fastq_bytes": "100;200",
                "analysis_include": "true",
            }
        )
        sample = blank(SAMPLE_COLUMNS)
        sample.update(
            {
                "study_id": study_id,
                "bioproject": bioproject,
                "sample_id": f"sample_{study_id}_{index}",
                "run_accession": accession,
                "fermentation_batch": f"{study_id}::batch_1",
                "fermentation_hours": str(hour),
                "fermentation_stage": stages[min(index, 2)],
                "analysis_include": "true",
            }
        )
        runs.append(run)
        samples.append(sample)
    return runs, samples


class PilotSelectionTests(unittest.TestCase):
    def test_more_complete_temporal_series_wins_deterministically(self) -> None:
        first = fixture_study("study_1", "PRJNA111111")
        second = fixture_study("study_2", "PRJNA222222")
        first_runs, first_samples = fixture_rows("study_1", "PRJNA111111", [0, 24, 48, 72])
        second_runs, second_samples = fixture_rows("study_2", "PRJNA222222", [0, 48, 72])

        rows = build_pilot_selection(
            [second, first], second_runs + first_runs, second_samples + first_samples, RULES
        )

        self.assertEqual(rows[0]["study_id"], "study_1")
        self.assertEqual(rows[0]["selected_as_pilot"], "true")
        self.assertEqual(rows[0]["selected_fastq_bytes"], "1200")
        self.assertEqual(rows[1]["selection_rank"], "2")

    def test_pending_study_is_compared_but_not_eligible(self) -> None:
        study = fixture_study("study_1", "PRJNA111111")
        study["include"] = "pending"
        runs, samples = fixture_rows("study_1", "PRJNA111111", [0, 24, 72])
        for row in runs + samples:
            row["analysis_include"] = "pending"

        rows = build_pilot_selection([study], runs, samples, RULES)

        self.assertEqual(rows[0]["eligible_primary"], "false")
        self.assertIn("study decision is pending", rows[0]["eligibility_reason"])
        self.assertEqual(rows[0]["candidate_runs"], "3")

    def test_malformed_fastq_bytes_fail_explicitly(self) -> None:
        study = fixture_study("study_1", "PRJNA111111")
        runs, samples = fixture_rows("study_1", "PRJNA111111", [0, 24, 72])
        runs[0]["fastq_bytes"] = "not-a-number"

        with self.assertRaisesRegex(ValueError, "Invalid FASTQ byte count"):
            build_pilot_selection([study], runs, samples, RULES)


if __name__ == "__main__":
    unittest.main()
