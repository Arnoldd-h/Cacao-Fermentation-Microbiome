"""Synthetic manifests verify inclusive scope; no network or sequencing data."""

import copy
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.config import load_json_yaml
from cacao_inventory.study_scope import build_study_manifest, build_stream_log, compare_live_runs, validate_scope


class StudyScopeTests(unittest.TestCase):
    def setUp(self):
        self.config = {"study_id": "synthetic_study", "bioproject": "PRJNA1", "expected_runs": 2, "expected_paired_runs": 2, "expected_batches": 1}
        self.runs = []
        self.samples = []
        for index in (1, 2):
            run = "SRR" + str(index)
            common = {"study_id": "synthetic_study", "bioproject": "PRJNA1", "run_accession": run, "analysis_include": "true", "sample_alias": "Bact" + str(index)}
            self.runs.append({**common, "library_layout": "PAIRED", "library_strategy": "AMPLICON", "sequencing_platform": "ILLUMINA",
                "fastq_ftp": ";".join("example.invalid/" + run + "_" + d + ".fastq.gz" for d in ("1", "2")),
                "fastq_md5": "a" * 32 + ";" + "b" * 32, "fastq_bytes": "10;20", "read_count": "2", "base_count": "100"})
            self.samples.append({**common, "sample_id": "sample" + str(index), "fermentation_batch": "batch1", "fermentation_hours": str((index - 1) * 12),
                "relative_time": str(index - 1), "fermentation_stage": "early" if index == 1 else "late", "sampling_stratum": "1",
                "marker": "16S rRNA", "country": "", "time_source": "synthetic fixture"})
        self.candidate = {"candidate_rules": [{"field": "sample_alias", "pattern": "^Bact"}, {"field": "library_strategy", "pattern": "^AMPLICON$"}]}

    def test_all_included_runs_preserve_known_and_unknown_metadata(self):
        rows = build_study_manifest(self.config, self.runs, self.samples)
        self.assertEqual([row["run_accession"] for row in rows], ["SRR1", "SRR2"])
        self.assertEqual(sum(int(row["estimated_bytes_total"]) for row in rows), 60)
        self.assertEqual(rows[0]["country"], "")
        self.assertEqual(rows[1]["fermentation_hours"], "12")
        self.assertEqual(rows[0]["time_source"], "synthetic fixture")

    def test_missing_duplicate_or_pending_samples_fail_without_subset(self):
        for samples in (self.samples[:1], [self.samples[0], self.samples[0]], [{**row, "analysis_include": "pending"} for row in self.samples]):
            with self.subTest(samples=samples), self.assertRaises(ValueError):
                build_study_manifest(self.config, self.runs, samples)

    def test_malformed_downloads_and_incompatible_technology_fail(self):
        for field, value in (("fastq_md5", "not-a-checksum"), ("fastq_bytes", "10"), ("fastq_ftp", ""), ("library_layout", "SINGLE"), ("sequencing_platform", "PACBIO_SMRT")):
            changed = copy.deepcopy(self.runs)
            changed[0][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                build_study_manifest(self.config, changed, self.samples)

    def test_changed_batch_count_and_duplicate_sample_ids_fail(self):
        for field, value in (("fermentation_batch", "different"), ("sample_id", "sample1"), ("sampling_stratum", "")):
            changed = copy.deepcopy(self.samples)
            changed[1][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                build_study_manifest(self.config, self.runs, changed)

    def test_live_sources_match_or_fail_on_new_missing_or_changed_runs(self):
        self.assertEqual(compare_live_runs(self.config, self.candidate, self.runs, self.runs)["eligible_runs"], 2)
        changed = copy.deepcopy(self.runs)
        changed[0]["fastq_md5"] = "c" * 32 + ";" + "b" * 32
        for live in (self.runs[:1], self.runs + [{**self.runs[0], "run_accession": "SRR3"}], changed, self.runs + self.runs):
            with self.subTest(live=live), self.assertRaises(ValueError):
                compare_live_runs(self.config, self.candidate, self.runs, live)

    def test_orphan_streams_are_logged_and_missing_mate_is_explicit(self):
        changed = copy.deepcopy(self.runs)
        changed[0]["fastq_ftp"] += ";example.invalid/SRR1.fastq.gz"
        changed[0]["fastq_md5"] += ";" + "c" * 32
        changed[0]["fastq_bytes"] += ";5"
        rows = build_study_manifest(self.config, changed, self.samples)
        self.assertEqual(rows[0]["estimated_bytes_total"], "30")
        self.assertIn("SRR1.fastq.gz", rows[0]["source_fastq_ftp"])
        log = build_stream_log(self.config, changed, self.samples)
        self.assertEqual([row["reason"] for row in log if row["decision"] == "exclude"], ["unpaired_stream"])
        changed[1].update(fastq_ftp="example.invalid/SRR2.fastq.gz", fastq_md5="b" * 32, fastq_bytes="20")
        with self.assertRaises(ValueError):
            build_study_manifest(self.config, changed, self.samples)
        cfg = {**self.config, "expected_paired_runs": 1}
        self.assertEqual(len(build_study_manifest(cfg, changed, self.samples)), 1)
        self.assertEqual(sum(row["reason"] == "no_complete_paired_files" for row in build_stream_log(cfg, changed, self.samples)), 1)

    def test_registered_paths_and_resource_budget_reject_unsafe_outputs(self):
        config = load_json_yaml(ROOT / "config/full_study.yaml")
        validate_scope(config, ROOT)
        for field, value in (("interim_dir", "data/raw/full"), ("dada2_dir", "results/dada2/pilot"), ("manifest", "../outside.tsv"), ("scope", "select_significant_runs")):
            changed = copy.deepcopy(config)
            changed[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                validate_scope(changed, ROOT)


if __name__ == "__main__":
    unittest.main()
