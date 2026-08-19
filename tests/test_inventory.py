from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_project_configuration  # noqa: E402
from cacao_inventory.inventory import build_samples  # noqa: E402


def run_row(accession: str, alias: str, **overrides: str) -> dict[str, str]:
    row = {
        "run_accession": accession,
        "sample_accession": f"SAMN{accession[3:]}",
        "sample_alias": alias,
        "sample_title": "fermentation bacteria",
        "sample_description": "",
        "library_strategy": "AMPLICON",
        "library_layout": "PAIRED",
        "instrument_platform": "ILLUMINA",
        "instrument_model": "Illumina MiSeq",
        "country": "Colombia:Antioquia",
    }
    row.update(overrides)
    return row


class InventoryTransformTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config, candidates = load_project_configuration(
            ROOT / "config" / "config.yaml", ROOT / "config" / "datasets.yaml"
        )
        cls.candidates = {candidate["bioproject"]: candidate for candidate in candidates}

    def test_group_max_duration_and_relative_time(self) -> None:
        rows = [
            run_row("SRR1000001", "BactAnI000z1"),
            run_row("SRR1000002", "BactAnI060z1"),
            run_row("SRR1000003", "BactAnI120z1"),
            run_row("SRR1000004", "YeastAnI120z1"),
        ]
        samples, selected = build_samples(
            self.candidates["PRJNA492720"], rows, self.config["temporal_stages"]
        )
        by_run = {row["run_accession"]: row for row in samples}
        self.assertEqual(selected, {"SRR1000001", "SRR1000002", "SRR1000003"})
        self.assertEqual(by_run["SRR1000001"]["relative_time"], "0")
        self.assertEqual(by_run["SRR1000002"]["relative_time"], "0.5")
        self.assertEqual(by_run["SRR1000003"]["relative_time"], "1")
        self.assertEqual(by_run["SRR1000002"]["fermentation_stage"], "mid")
        self.assertEqual(by_run["SRR1000003"]["fermentation_duration_hours"], "120")

    def test_unknown_duration_does_not_invent_relative_time(self) -> None:
        rows = [
            run_row(
                "SRR2000001",
                "04.4.5.1day.16S.b",
                country="USA: Hawaii, Oahu",
            )
        ]
        samples, _ = build_samples(
            self.candidates["PRJNA865318"], rows, self.config["temporal_stages"]
        )
        self.assertEqual(samples[0]["fermentation_hours"], "24")
        self.assertEqual(samples[0]["fermentation_duration_hours"], "")
        self.assertEqual(samples[0]["relative_time"], "")
        self.assertEqual(samples[0]["fermentation_stage"], "")

    def test_non_inoculated_invitro_time_is_parsed(self) -> None:
        rows = [
            run_row(
                "SRR3000001",
                "invitro16S_01",
                sample_description="metagenome of non-inoculated fermentation mass at 48hr-R3",
                country="Colombia",
            )
        ]
        samples, _ = build_samples(
            self.candidates["PRJNA1104253"], rows, self.config["temporal_stages"]
        )
        self.assertEqual(samples[0]["fermentation_hours"], "48")
        self.assertEqual(samples[0]["relative_time"], "0.5")
        self.assertEqual(samples[0]["replicate"], "3")
        self.assertEqual(samples[0]["analysis_include"], "pending")


if __name__ == "__main__":
    unittest.main()
