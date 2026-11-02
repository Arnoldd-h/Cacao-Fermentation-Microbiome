"""Synthetic cohort fixtures exercise DADA2 membership and promotion guards."""

from copy import deepcopy
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("study_dada2_cohort", ROOT / "scripts/study_dada2/cohort.py")
COHORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COHORT)


def manifest_rows():
    return [{
        "study_id": "synthetic_study", "sample_id": f"sample_{index}", "run_accession": f"run_{index}",
        "pilot_order": str(index), "fermentation_batch": batch, "fermentation_hours": hours,
        "relative_time": relative, "fermentation_duration_hours": "24", "estimated_bytes_total": "1000",
        "sampling_stratum": "", "replicate": "01", "notes": "Synthetic test data",
    } for index, (batch, hours, relative) in enumerate((("batch_a", "0", "0"), ("batch_a", "24", "1"), ("batch_b", "0", "0")), 1)]


def terminal_rows(manifest):
    return [{**{field: sample[field] for field in COHORT.IDENTITY_FIELDS}, "direction": direction, "pairs": str(index * 100)}
            for index, sample in enumerate(manifest, 1) for direction in ("R1", "R2")]


class QualityRowsTests(unittest.TestCase):
    def setUp(self):
        self.manifest = manifest_rows()
        self.terminal = terminal_rows(self.manifest)

    def test_output_order_and_exact_identities_follow_manifest(self):
        rows = COHORT.quality_rows(self.manifest, list(reversed(self.terminal)))
        self.assertEqual([(row["sample_id"], row["read_direction"]) for row in rows],
                         [(sample["sample_id"], direction) for sample in self.manifest for direction in ("R1", "R2")])
        self.assertEqual(rows[0], {"study_id": "synthetic_study", "sample_id": "sample_1", "run_accession": "run_1", "read_direction": "R1", "total_sequences": 100})

    def test_duplicate_missing_extra_or_invalid_mate_fails(self):
        cases = [self.terminal + [self.terminal[0]], self.terminal[:-1],
                 self.terminal + [{**self.terminal[0], "sample_id": "unexpected"}],
                 [{**self.terminal[0], "direction": "F"}, *self.terminal[1:]]]
        for rows in cases:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                COHORT.quality_rows(self.manifest, rows)

    def test_duplicate_manifest_sample_or_run_fails(self):
        for field in ("sample_id", "run_accession"):
            manifest = deepcopy(self.manifest)
            manifest[1][field] = manifest[0][field]
            with self.subTest(field=field), self.assertRaises(ValueError):
                COHORT.quality_rows(manifest, self.terminal)

    def test_identity_count_mismatch_or_invalid_count_fails(self):
        for field, value in (("study_id", "other"), ("run_accession", "other"),
                             ("pairs", "99"), ("pairs", "1.5"), ("pairs", "-1"), ("pairs", True)):
            rows = deepcopy(self.terminal)
            rows[0][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                COHORT.quality_rows(self.manifest, rows)

    def test_zero_counts_are_preserved_for_pre_execution_quality_rows(self):
        for row in self.terminal:
            row["pairs"] = "0"
        self.assertEqual([row["total_sequences"] for row in COHORT.quality_rows(self.manifest, self.terminal)], [0] * 6)


class ValidateCohortTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.manifest = manifest_rows()
        self.terminal = terminal_rows(self.manifest)
        self.summary = {"study_id": "synthetic_study", "runs": 3, "pairs": 600, "pairs_below_min_len": 10}
        self.policy = {"expected_runs": 3, "expected_batches": 2, "expected_batch_times": 3}
        self.tables = {"sample_metadata.tsv": deepcopy(self.manifest), "read_tracking.tsv": [], "merge_diagnostics.tsv": []}
        for index, sample in enumerate(self.manifest, 1):
            identities = {field: sample[field] for field in COHORT.IDENTITY_FIELDS}
            self.tables["read_tracking.tsv"].append({
                **identities, "fermentation_hours": sample["fermentation_hours"],
                "input": str(100 * index), "raw_input": str(100 * index), "filtered": str(80 * index),
                "denoised_f": str(75 * index), "denoised_r": str(74 * index),
                "merged": str(70 * index), "nonchim": str(60 * index),
            })
            self.tables["merge_diagnostics.tsv"].append({
                "study_id": sample["study_id"], "sample_id": sample["sample_id"],
                "accepted_pairs": str(70 * index), "rejected_pairs": str(4 * index),
            })
        self.tables["error_learning.tsv"] = [{
            "read_direction": direction, "self_consistency_rounds": "3", "maximum_rounds": "10",
            "repeated_error_matrix": "TRUE", "reached_iteration_limit": "FALSE", "converged_before_limit": "TRUE",
            "transition_base_observations": "1000",
        } for direction in ("F", "R")]
        self.write_tables()

    def write_tables(self):
        for name, rows in self.tables.items():
            with (self.directory / name).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

    def validate(self):
        return COHORT.validate_cohort(self.directory, self.manifest, self.terminal, self.summary, self.policy)

    def test_complete_cohort_ready_for_taxonomy_and_counts_reconcile(self):
        result = self.validate()
        self.assertEqual(result["status"], "valid")
        self.assertEqual(result["promotion_status"], "ready_for_taxonomy")
        self.assertEqual((result["samples"], result["batches"], result["batch_times"]), (3, 2, 3))
        self.assertEqual(result["unknown_sampling_stratum_samples"], 3)
        self.assertEqual((result["input_pairs"], result["filtered_pairs"], result["merged_pairs"], result["nonchim_pairs"]), (600, 480, 420, 360))
        self.assertEqual(result["pairs_at_least_min_len"], 590)
        self.assertEqual(result["minimum_nonchim_pairs_per_sample"], 60)
        self.assertEqual(result["retention"]["nonchim_fraction_of_input"], 0.6)

    def test_numeric_formatting_changes_only_are_accepted(self):
        metadata = self.tables["sample_metadata.tsv"][0]
        for field in COHORT.NUMERIC_METADATA:
            metadata[field] = metadata[field] + ".0"
        self.tables["read_tracking.tsv"][0]["fermentation_hours"] = "0.0"
        self.write_tables()
        self.assertEqual(self.validate()["promotion_status"], "ready_for_taxonomy")

    def test_metadata_alteration_including_unknown_strata_fails(self):
        for field, value in (("sampling_stratum", "1"), ("notes", "Altered"), ("fermentation_hours", "0.01"), ("replicate", "1")):
            with self.subTest(field=field):
                self.tables["sample_metadata.tsv"] = deepcopy(self.manifest)
                self.tables["sample_metadata.tsv"][0][field] = value
                self.write_tables()
                with self.assertRaisesRegex(ValueError, "metadata altered"):
                    self.validate()

    def test_missing_or_added_manifest_columns_fail(self):
        for alteration in ("missing", "added"):
            self.tables["sample_metadata.tsv"] = deepcopy(self.manifest)
            for row in self.tables["sample_metadata.tsv"]:
                if alteration == "missing":
                    del row["notes"]
                else:
                    row["invented_field"] = ""
            self.write_tables()
            with self.subTest(alteration=alteration), self.assertRaisesRegex(ValueError, "every manifest field"):
                self.validate()

    def test_duplicate_missing_or_extra_cohort_members_fail(self):
        original = deepcopy(self.tables)
        for name in ("sample_metadata.tsv", "read_tracking.tsv", "merge_diagnostics.tsv"):
            for kind in ("duplicate", "missing", "extra"):
                self.tables = deepcopy(original)
                rows = self.tables[name]
                if kind == "duplicate":
                    rows.append(deepcopy(rows[0]))
                elif kind == "missing":
                    rows.pop()
                else:
                    rows.append({**rows[0], "sample_id": "unexpected"})
                self.write_tables()
                with self.subTest(name=name, kind=kind), self.assertRaises(ValueError):
                    self.validate()

    def test_tracking_identity_input_or_zero_nonchim_fails(self):
        original = deepcopy(self.tables)
        for field, value in (("run_accession", "wrong"), ("study_id", "wrong"), ("input", "101"),
                             ("raw_input", "101"), ("nonchim", "0"), ("denoised_f", "999")):
            self.tables = deepcopy(original)
            self.tables["read_tracking.tsv"][0][field] = value
            self.write_tables()
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.validate()

    def test_merge_counts_or_study_mismatch_fails(self):
        for field, value in (("accepted_pairs", "71"), ("study_id", "wrong"), ("rejected_pairs", "-1")):
            original = deepcopy(self.tables["merge_diagnostics.tsv"])
            self.tables["merge_diagnostics.tsv"][0][field] = value
            self.write_tables()
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.validate()
            self.tables["merge_diagnostics.tsv"] = original

    def test_filtered_counts_cannot_exceed_length_eligible_pairs(self):
        self.summary["pairs_below_min_len"] = 121
        with self.assertRaisesRegex(ValueError, "minimum-length eligibility"):
            self.validate()

    def test_terminal_summary_count_or_identity_mismatch_fails(self):
        original = deepcopy(self.summary)
        for field, value in (("study_id", "wrong"), ("runs", 2), ("pairs", 601), ("pairs_below_min_len", 601)):
            self.summary = {**original, field: value}
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.validate()

    def test_expected_design_counts_are_enforced(self):
        for field in ("expected_runs", "expected_batches", "expected_batch_times"):
            original = self.policy[field]
            self.policy[field] = original + 1
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, field):
                self.validate()
            self.policy[field] = original

    def test_nonconverged_valid_outputs_need_review_without_throwing(self):
        self.tables["error_learning.tsv"][1].update({
            "self_consistency_rounds": "10", "repeated_error_matrix": "FALSE",
            "reached_iteration_limit": "TRUE", "converged_before_limit": "FALSE",
        })
        self.write_tables()
        result = self.validate()
        self.assertEqual(result["status"], "valid")
        self.assertEqual(result["promotion_status"], "needs_review")
        self.assertIn("R", result["promotion_reason"])
        self.assertFalse(result["error_learning"]["R"]["converged_before_limit"])

    def test_missing_error_direction_or_inconsistent_rounds_flags_fail(self):
        original = deepcopy(self.tables["error_learning.tsv"])
        mutations = (("read_direction", "F"), ("self_consistency_rounds", "0"),
                     ("self_consistency_rounds", "11"), ("maximum_rounds", "0"),
                     ("reached_iteration_limit", "TRUE"), ("converged_before_limit", "FALSE"),
                     ("repeated_error_matrix", "FALSE"), ("converged_before_limit", "true"))
        for field, value in mutations:
            self.tables["error_learning.tsv"] = deepcopy(original)
            self.tables["error_learning.tsv"][1][field] = value
            self.write_tables()
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                self.validate()
        self.tables["error_learning.tsv"] = original[:1]
        self.write_tables()
        with self.assertRaisesRegex(ValueError, "exactly F and R"):
            self.validate()


if __name__ == "__main__":
    unittest.main()
