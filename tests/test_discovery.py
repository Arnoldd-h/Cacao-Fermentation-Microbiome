from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.discovery import preliminary_screen  # noqa: E402


class DiscoveryTriageTests(unittest.TestCase):
    def test_wgs_only_is_outside_phase_one(self) -> None:
        summary = {"project_title": "Cocoa fermentation shotgun"}
        runs = [{"library_strategy": "WGS", "sample_title": "fermentation"}]
        screen, _, _, _ = preliminary_screen(summary, runs, configured=False)
        self.assertEqual(screen, "exclude_phase_i_wgs")

    def test_temporal_16s_amplicon_is_manual_priority(self) -> None:
        summary = {"project_title": "Cocoa fermentation 16S succession"}
        runs = [
            {
                "library_strategy": "AMPLICON",
                "sample_title": "bacterial fermentation at 48 hr",
                "sample_alias": "sample_48hr_16S",
            }
        ]
        screen, _, bacterial, temporal = preliminary_screen(summary, runs, configured=False)
        self.assertEqual(screen, "manual_review_priority")
        self.assertTrue(bacterial)
        self.assertTrue(temporal)

    def test_configured_projects_are_not_reclassified(self) -> None:
        summary = {"project_title": "Cocoa fermentation"}
        screen, _, _, _ = preliminary_screen(summary, [], configured=True)
        self.assertEqual(screen, "configured")


if __name__ == "__main__":
    unittest.main()
