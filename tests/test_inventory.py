from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_project_configuration  # noqa: E402
from cacao_inventory.inventory import build_runs, build_samples  # noqa: E402


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
        candidate = {
            "study_id": "fixture_unknown_duration", "bioproject": "PRJNA999999",
            "screening_status": "pending", "duration_mode": "unknown",
            "candidate_rules": [{"field": "sample_alias", "pattern": "16S"}],
            "time_parser": {"field": "sample_alias", "pattern": r"(?P<days>\d+)day", "unit": "days"},
        }
        samples, _ = build_samples(candidate, rows, self.config["temporal_stages"])
        self.assertEqual(samples[0]["fermentation_hours"], "24")
        self.assertEqual(samples[0]["fermentation_duration_hours"], "")
        self.assertEqual(samples[0]["relative_time"], "")
        self.assertEqual(samples[0]["fermentation_stage"], "")
        self.assertEqual(samples[0]["fermentation_batch"], "")

    def test_resolved_invitro_project_does_not_enter_primary_candidates(self) -> None:
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
        self.assertEqual(samples, [])
        self.assertEqual(self.candidates["PRJNA1104253"]["screening_status"], "exclude")

    def test_costa_rica_2017_selector_requires_bacterial_submitted_files(self) -> None:
        rows = [
            run_row("ERR9000001", "F1T20", submitted_ftp="host/F1T20_F_16S.fastq.gz;host/F1T20_R_16S.fastq.gz"),
            run_row("ERR9000002", "F1T20", submitted_ftp="host/F1T20_F_ITS.fastq.gz;host/F1T20_R_ITS.fastq.gz"),
            run_row("ERR9000003", "F3T20", submitted_ftp="host/F3T20_F_16S.fastq.gz;host/F3T20_R_16S.fastq.gz"),
            run_row("ERR9000004", "F1T20"),
        ]
        samples, selected = build_samples(self.candidates["PRJEB40850"], rows, self.config["temporal_stages"])
        self.assertEqual(selected, {"ERR9000001"})
        self.assertEqual(samples[0]["fermentation_hours"], "20")
        self.assertEqual(samples[0]["fermentation_duration_hours"], "92")
        self.assertTrue(samples[0]["fermentation_batch"].endswith("::1"))

    def test_costa_rica_2019_selector_separates_pacbio_from_its_and_wgs(self) -> None:
        rows = [
            run_row("ERR9000011", "F01T120", instrument_model="Sequel II", instrument_platform="PACBIO_SMRT", library_layout="SINGLE", submitted_ftp="host/F01T120.16S.fastq.gz"),
            run_row("ERR9000012", "F01T120", submitted_ftp="host/F01T120.ITS.fastq.gz"),
            run_row("ERR9000013", "F01T120", library_strategy="WGS", instrument_model="Illumina NovaSeq 6000"),
        ]
        samples, selected = build_samples(self.candidates["PRJEB57747"], rows, self.config["temporal_stages"])
        self.assertEqual(selected, {"ERR9000011"})
        self.assertEqual(samples[0]["relative_time"], "1")
        self.assertEqual(samples[0]["region_16s"], "full-length")

    def test_conflicting_cameroon_144_hour_record_is_preserved_pending(self) -> None:
        row = run_row("SRR9000021", "J_conflict", sample_title="Bacteria_metagenome_from_cocoa_beans_fermentation_time_144_fermentation_type_Heap_Control_second_replicate")
        samples, selected = build_samples(self.candidates["PRJNA420946"], [row], self.config["temporal_stages"])
        self.assertEqual(selected, {"SRR9000021"})
        self.assertEqual(samples[0]["fermentation_hours"], "144")
        self.assertEqual(samples[0]["analysis_include"], "pending")
        self.assertEqual(samples[0]["relative_time"], "")

    def test_submitted_files_are_preserved_as_observed_run_metadata(self) -> None:
        row = run_row("ERR9000031", "F1T0", submitted_ftp="host/F1T0_F_16S.fastq.gz", submitted_format="FASTQ")
        runs, _ = build_runs(self.candidates["PRJEB40850"], [row], {"ERR9000031"})
        self.assertEqual(runs[0]["submitted_ftp"], row["submitted_ftp"])
        self.assertEqual(runs[0]["submitted_format"], "FASTQ")

    def test_raw_label_maps_to_zero_and_preserves_sample_variety(self) -> None:
        candidate = {
            "bioproject": "PRJNA999999",
            "study_id": "fixture_mexico",
            "screening_status": "include",
            "country": "Mexico",
            "region": "Tabasco",
            "cacao_variety": "multiple",
            "region_16s": "V3-V4",
            "season": "",
            "candidate_rules": [
                {"field": "sample_alias", "pattern": "^16S_"},
                {
                    "field": "sample_title",
                    "pattern": "(?P<variety>Criollo|Forastero).*Batch_(?P<batch>\\d+).*replicate_(?P<replicate>\\d+)",
                },
            ],
            "time_parser": {
                "field": "sample_title",
                "pattern": "(?:(?P<time_zero>raw)|after (?P<hours>\\d+)h of fermentation)",
                "unit": "hours",
                "zero_group": "time_zero",
            },
            "duration_mode": "fixed",
            "duration_hours": 120,
            "duration_group_fields": ["variety", "batch"],
            "value_maps": {},
        }
        rows = [
            run_row(
                "SRR4000001",
                "16S_fixture_1",
                sample_title=(
                    "bacterial 16S under traditional conditions of raw Criollo cocoa beans "
                    "from Mexico Batch_1_ replicate_2"
                ),
            )
        ]
        samples, _ = build_samples(candidate, rows, self.config["temporal_stages"])
        self.assertEqual(samples[0]["fermentation_hours"], "0")
        self.assertEqual(samples[0]["relative_time"], "0")
        self.assertEqual(samples[0]["fermentation_stage"], "early")
        self.assertEqual(samples[0]["cacao_variety"], "Criollo")
        self.assertEqual(samples[0]["replicate"], "2")

    def test_spontaneous_comparator_is_not_labeled_as_technical_control(self) -> None:
        candidate = {
            "bioproject": "PRJEB99999",
            "study_id": "fixture_ecuador",
            "screening_status": "pending",
            "screening_reason": "Primer verification pending",
        }
        row = run_row(
            "ERR4000001",
            "F1T0",
            sample_description="Fermentation 1, 0 h, negative control",
        )
        runs, exclusions = build_runs(candidate, [row], set())
        self.assertEqual(runs[0]["marker_classification"], "not_verified_16s_subset")
        self.assertEqual(exclusions[0]["metric"], "not_verified_16s_subset")

    def test_mixed_16s_its_description_stays_unresolved(self) -> None:
        candidate = {
            "bioproject": "PRJEB99998",
            "study_id": "fixture_mixed_marker",
            "screening_status": "pending",
            "screening_reason": "Run-level marker mapping pending",
        }
        row = run_row(
            "ERR4000002",
            "mixed_marker",
            sample_description="V4 region of bacterial 16S rRNA gene and fungal ITS1 region",
        )
        runs, _ = build_runs(candidate, [row], set())
        self.assertEqual(runs[0]["marker_classification"], "not_verified_16s_subset")


if __name__ == "__main__":
    unittest.main()
