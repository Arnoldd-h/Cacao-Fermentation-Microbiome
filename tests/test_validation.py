from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_project_configuration  # noqa: E402
from cacao_inventory.schema import (  # noqa: E402
    EXCLUSION_COLUMNS,
    RUN_COLUMNS,
    SAMPLE_COLUMNS,
    STUDY_COLUMNS,
)
from cacao_inventory.validation import InventoryValidationError, validate_inventory  # noqa: E402


def blank(columns: list[str]) -> dict[str, str]:
    return {column: "" for column in columns}


def valid_tables() -> dict[str, list[dict[str, str]]]:
    study = blank(STUDY_COLUMNS)
    study.update(
        {
            "study_id": "fixture_study",
            "bioproject": "PRJNA123456",
            "marker": "16S rRNA",
            "include": "true",
            "raw_data_available": "true",
            "number_samples": "3",
        }
    )
    runs = []
    samples = []
    for index, (hours, relative, stage) in enumerate(
        [("0", "0", "early"), ("48", "0.5", "mid"), ("96", "1", "late")], start=1
    ):
        run_accession = f"SRR900000{index}"
        run = blank(RUN_COLUMNS)
        run.update(
            {
                "study_id": "fixture_study",
                "bioproject": "PRJNA123456",
                "run_accession": run_accession,
                "analysis_include": "true",
            }
        )
        sample = blank(SAMPLE_COLUMNS)
        sample.update(
            {
                "study_id": "fixture_study",
                "sample_id": f"fixture_{index}",
                "run_accession": run_accession,
                "bioproject": "PRJNA123456",
                "fermentation_hours": hours,
                "fermentation_duration_hours": "96",
                "relative_time": relative,
                "fermentation_stage": stage,
                "analysis_include": "true",
            }
        )
        runs.append(run)
        samples.append(sample)
    return {"studies": [study], "runs": runs, "samples": samples, "exclusion_log": []}


class InventoryValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config, _ = load_project_configuration(
            ROOT / "config" / "config.yaml", ROOT / "config" / "datasets.yaml"
        )

    def test_valid_inventory_passes(self) -> None:
        summary = validate_inventory(valid_tables(), self.config)
        self.assertEqual(summary["samples"], 3)

    def test_incorrect_relative_time_fails(self) -> None:
        tables = valid_tables()
        tables["samples"][1]["relative_time"] = "0.4"
        with self.assertRaisesRegex(InventoryValidationError, "Incorrect relative_time"):
            validate_inventory(tables, self.config)

    def test_excluded_run_requires_log_entry(self) -> None:
        tables = valid_tables()
        tables["runs"][0]["analysis_include"] = "false"
        with self.assertRaisesRegex(InventoryValidationError, "missing from exclusion_log"):
            validate_inventory(tables, self.config)

        exclusion = blank(EXCLUSION_COLUMNS)
        exclusion.update(
            {
                "entity_type": "run",
                "entity_id": tables["runs"][0]["run_accession"],
                "decision": "exclude",
            }
        )
        tables["exclusion_log"].append(exclusion)
        validate_inventory(tables, self.config)


if __name__ == "__main__":
    unittest.main()
