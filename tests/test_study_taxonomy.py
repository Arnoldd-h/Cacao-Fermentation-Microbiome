"""Synthetic study-taxonomy gates; classification subprocesses are never run."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("independent_study_taxonomy_wrapper", ROOT / "scripts/study_taxonomy/execute.py")
TAX = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TAX)


class StudyTaxonomyTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        (self.root / "config").mkdir()
        self.profile = json.loads((ROOT / "config/taxonomy.yaml").read_text())
        self.profile_path = self.root / "config/taxonomy.yaml"
        self.profile_path.write_text(json.dumps(self.profile))
        self.policy = {
            "schema_version": "1.0.0", "method": "independent_study_registered_taxonomy",
            "source_dada2_policy": "config/dada2.json", "source_primer_review_policy": "config/primers.json",
            "taxonomy_profile": "config/taxonomy.yaml", "taxonomy_profile_sha256": TAX.sha256_file(self.profile_path),
            "output_dir": "results/taxonomy/synthetic/primary", "audit_dir": "results/taxonomy/synthetic/audit",
            "repeat_dir": "results/intermediate/taxonomy/synthetic/repeat", "repeat_threads": 1,
            "require_identical_repetition": True, "expected_samples": 2, "expected_asvs": 3, "expected_nonchim_pairs": 10,
            "threads": 2, "memory_mb": 4500, "require_dada2_promotion": "ready_for_taxonomy",
            "residual_primer_policy": "retain_all_diagnosed_asvs", "sample_exclusions": False,
            "asv_exclusions": False, "parameter_overrides": False, "processing_ceiling": "taxonomy",
        }
        self.policy_path = self.root / "config/study_taxonomy.json"
        self.write_policy()
        self.review = {"source_dada2_policy": "config/dada2.json", "output_dir": "results/diagnostics/synthetic/primers"}
        self.dada = {"output_dir": "results/dada2/synthetic/run", "audit_dir": "results/dada2/synthetic/audit",
                     "expected_runs": 2, "expected_dada2_version": "1.34.0", "expected_r_version": "4.4.3",
                     "source_terminal_policy": "config/terminal.json"}
        self.scope = {"study_id": "synthetic", "expected_runs": 2, "taxonomy_dir": "results/taxonomy/synthetic/full",
                      "resources": {"memory_budget_mb": 4500}, "manifest": "metadata/manifest.tsv"}
        self.processing = {"dada2": {"random_seed": self.profile["classification"]["random_seed"]}}
        (self.root / "config/terminal.json").write_text(json.dumps({"output_dir": "results/qc/synthetic/terminal"}))
        self.source_paths = [self.root / "config/dada2.json"]
        self.loader = patch.object(TAX.REVIEW, "load_policy", return_value=(self.review, self.dada, self.scope, self.processing, self.source_paths))
        self.loader.start()
        self.addCleanup(self.loader.stop)
        self.native = {"status": "valid", "source_git_commit": "a" * 40, "study_id": "synthetic", "samples": 2,
                       "asvs": 3, "total_reads": 10, "excluded_asvs": 0}

    def write_policy(self):
        self.policy_path.write_text(json.dumps(self.policy))

    def context(self):
        return TAX.load_policy(self.policy_path, self.root)

    def test_policy_retains_immutable_profile_and_independent_paths(self):
        context = self.context()
        self.assertEqual(context[0], self.policy)
        self.assertEqual(context[1], self.profile)
        self.assertEqual(context[4]["study_id"], "synthetic")

    def test_profile_bytes_seed_or_missing_native_adapter_fail(self):
        self.profile_path.write_text(json.dumps(self.profile) + "\n")
        with self.assertRaisesRegex(ValueError, "profile hash"):
            self.context()
        for field in ("random_seed", "native_rng"):
            profile = deepcopy(self.profile)
            if field == "random_seed": profile["classification"][field] += 1
            else: profile["classification"].pop(field)
            self.profile_path.write_text(json.dumps(profile))
            self.policy["taxonomy_profile_sha256"] = TAX.sha256_file(self.profile_path)
            self.write_policy()
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "adapter/version/seed"):
                self.context()

    def test_invalid_policy_flags_resources_or_source_binding_fail(self):
        original = deepcopy(self.policy)
        cases = (("asv_exclusions", True), ("parameter_overrides", True), ("require_identical_repetition", False),
                 ("expected_samples", 1), ("expected_asvs", True), ("memory_mb", 4501), ("threads", 1),
                 ("repeat_threads", 2), ("source_dada2_policy", "config/another.json"))
        for field, value in cases:
            self.policy = {**original, field: value}
            self.write_policy()
            with self.subTest(field=field), self.assertRaises(ValueError): self.context()

    def test_escape_overlap_baseline_and_wrong_repeat_namespace_fail(self):
        original = deepcopy(self.policy)
        for field, value in (("output_dir", "../outside"), ("output_dir", self.scope["taxonomy_dir"]),
                             ("audit_dir", self.policy["output_dir"]), ("audit_dir", self.policy["output_dir"] + "/audit"),
                             ("repeat_dir", "results/taxonomy/synthetic/repeat")):
            self.policy = {**original, field: value}
            self.write_policy()
            with self.subTest(field=field, value=value), self.assertRaises(ValueError): self.context()

    def test_inputs_list_binds_validators_and_workflow_without_hashing_reference(self):
        directory = self.root / self.dada["audit_dir"]
        directory.mkdir(parents=True)
        (directory / "provenance.json").write_text(json.dumps({"outputs": []}))
        with patch.object(TAX, "file_records", side_effect=AssertionError("No expensive hash during workflow parsing")):
            paths = TAX.inputs(self.policy_path, self.root)
        for relative in ("scripts/taxonomy/validate_taxonomy.py", "scripts/taxonomy/audit_repetition.py"):
            # Native source listing uses files present in the repository.
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# synthetic source\n")
        paths = TAX.inputs(self.policy_path, self.root)
        for relative in ("scripts/taxonomy/validate_taxonomy.py", "scripts/taxonomy/audit_repetition.py", "workflow/study_taxonomy.Snakefile", "scripts/study_taxonomy/execute.py"):
            self.assertIn(self.root / relative, paths)

    def test_unready_dada2_or_changed_cohort_counts_fail_before_classification(self):
        good = {"promotion_status": "ready_for_taxonomy", "native_validation": {"samples": 2, "asvs": 3, "nonchim_pairs": 10}}
        for field in ("promotion_status", "samples", "asvs", "nonchim_pairs"):
            source = deepcopy(good)
            if field == "promotion_status": source[field] = "needs_review"
            else: source["native_validation"][field] += 1
            with self.subTest(field=field), patch.object(TAX.REVIEW, "checked_source", return_value=(source, [])), self.assertRaises(ValueError):
                TAX.check_sources(self.context(), self.root)

    def setup_native(self, repeat=False):
        directory = self.root / self.policy["repeat_dir" if repeat else "output_dir"]
        directory.mkdir(parents=True, exist_ok=True)
        threads = self.policy["repeat_threads" if repeat else "threads"]
        source = self.root / self.dada["output_dir"]
        source.mkdir(parents=True, exist_ok=True)
        upstream = {"status": "valid", "source_git_commit": "b" * 40}
        (source / "validation.json").write_text(json.dumps(upstream))
        expected_paths = ["config/taxonomy.yaml", "config/config.yaml", self.profile["reference"]["path"], self.profile["reference"]["provenance_path"],
                          *[self.dada["output_dir"] + "/" + name for name in ("SUCCESS", "validation.json", "asv_sequences.tsv", "asv_counts.tsv", "sample_metadata.tsv")]]
        provenance = {"git_commit": "a" * 40, "effective_threads": threads, "classification": self.profile["classification"],
                      "reference": self.profile["reference"], "rng_kind": self.profile["classification"]["rng_kind"], "upstream_validation": upstream,
                      "inputs": [{"path": path} for path in expected_paths],
                      "command": ["Rscript", "scripts/taxonomy/assign_taxonomy.R", "--config", str(self.profile_path), "--input-dir", str(source),
                                  "--output-dir", str(directory), "--threads", str(threads)],
                      "started_at_utc": "2026-01-01T00:02:00+00:00" if repeat else "2026-01-01T00:00:00+00:00",
                      "finished_at_utc": "2026-01-01T00:03:00+00:00" if repeat else "2026-01-01T00:01:00+00:00"}
        (directory / "config_snapshot.yaml").write_text(json.dumps(self.profile))
        (directory / "provenance.json").write_text(json.dumps(provenance))
        (directory / "SUCCESS").write_text("a" * 40 + "\n")
        (directory / "software_versions.tsv").write_text("software\tversion\nR\t4.4.3\ndada2\t1.34.0\ntaxonomy_native_rng\tr_seeded_query_mt19937_v1\n")
        return directory, provenance

    def test_native_results_require_correct_producer_thread_source_and_versions(self):
        directory, original = self.setup_native()
        with patch.object(TAX.NATIVE, "validate", return_value=self.native):
            self.assertEqual(TAX.native_result(self.context(), self.root, expected_revision="a" * 40), self.native)
            for field, value in (("effective_threads", 1), ("git_commit", "c" * 40), ("rng_kind", ["wrong"]),
                                 ("upstream_validation", {"status": "wrong"})):
                (directory / "provenance.json").write_text(json.dumps({**original, field: value}))
                with self.subTest(field=field), self.assertRaises(ValueError):
                    TAX.native_result(self.context(), self.root, expected_revision="a" * 40)
            (directory / "provenance.json").write_text(json.dumps(original))
            (directory / "software_versions.tsv").write_text("software\tversion\nR\twrong\ndada2\t1.34.0\n")
            with self.assertRaisesRegex(ValueError, "software/seeded adapter"):
                TAX.native_result(self.context(), self.root)

    def test_native_command_and_input_manifest_cannot_substitute_another_cohort(self):
        directory, original = self.setup_native()
        with patch.object(TAX.NATIVE, "validate", return_value=self.native):
            for flag in ("--config", "--input-dir", "--output-dir", "--threads"):
                changed = deepcopy(original)
                changed["command"][changed["command"].index(flag) + 1] = "wrong"
                (directory / "provenance.json").write_text(json.dumps(changed))
                with self.subTest(flag=flag), self.assertRaisesRegex(ValueError, "command differs"):
                    TAX.native_result(self.context(), self.root)
            for mutation in ("missing", "substituted"):
                changed = deepcopy(original)
                item = next(row for row in changed["inputs"] if row["path"].endswith("sample_metadata.tsv"))
                if mutation == "missing": changed["inputs"].remove(item)
                else: item["path"] = "results/dada2/other/run/sample_metadata.tsv"
                (directory / "provenance.json").write_text(json.dumps(changed))
                with self.subTest(mutation=mutation), self.assertRaisesRegex(ValueError, "inputs do not bind"):
                    TAX.native_result(self.context(), self.root)

    def setup_repetition(self):
        baseline, _ = self.setup_native()
        repeat, _ = self.setup_native(repeat=True)
        for directory in (baseline, repeat):
            for name in TAX.REPETITION.TABLES:
                (directory / name).write_text("study_id\tasv_id\tvalue\nsynthetic\tA\t5\n")
        audit = self.root / self.policy["audit_dir"]
        audit.mkdir(parents=True)
        self.refresh_repetition_record()
        return baseline, repeat, audit

    def refresh_repetition_record(self):
        baseline, repeat = (self.root / self.policy[key] for key in ("output_dir", "repeat_dir"))
        executions = {}
        for label, directory in (("baseline", baseline), ("repeat", repeat)):
            provenance = json.loads((directory / "provenance.json").read_text())
            executions[label] = {"directory": directory.relative_to(self.root).as_posix(), "git_commit": "a" * 40,
                                 "config": self.profile, "provenance_sha256": TAX.sha256_file(directory / "provenance.json"),
                                 **{key: provenance[key] for key in ("effective_threads", "started_at_utc", "finished_at_utc")}}
        record = {"all_identical": True, "tables": TAX.REPETITION.compare_tables(baseline, repeat), "executions": executions}
        (self.root / self.policy["audit_dir"] / "repetition.json").write_text(json.dumps(record))

    def test_repetition_checks_identical_tables_and_chronology(self):
        _, repeat, _ = self.setup_repetition()
        self.assertEqual(TAX.check_repetition(self.context(), self.root)["tables_checked"], 7)
        (repeat / "taxonomy.tsv").write_text("study_id\tasv_id\tvalue\nsynthetic\tA\t6\n")
        with self.assertRaisesRegex(ValueError, "repetition differs"):
            TAX.check_repetition(self.context(), self.root)

    def test_repetition_input_drift_or_reversed_chronology_fails(self):
        _, repeat, _ = self.setup_repetition()
        original = json.loads((repeat / "provenance.json").read_text())
        for kind in ("inputs", "chronology"):
            changed = deepcopy(original)
            if kind == "inputs": changed["inputs"].append({"path": "different"})
            else: changed["started_at_utc"] = "2025-12-31T23:59:00+00:00"
            (repeat / "provenance.json").write_text(json.dumps(changed))
            self.refresh_repetition_record()
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                TAX.check_repetition(self.context(), self.root)

    def test_complete_and_partial_directories_are_never_overwritten(self):
        for key in ("output_dir", "audit_dir", "repeat_dir"):
            directory = self.root / self.policy[key]
            directory.mkdir(parents=True)
            marker = directory / "existing.txt"
            marker.write_text("preserve")
            with self.subTest(key=key), patch.object(TAX, "check_sources") as checker, self.assertRaisesRegex(ValueError, "must be preserved"):
                TAX.execute(self.policy_path, self.root)
            checker.assert_not_called()
            self.assertEqual(marker.read_text(), "preserve")
            marker.unlink()

    def test_execution_is_sequential_and_success_requires_final_validation(self):
        events = []
        audit = self.root / self.policy["audit_dir"]
        def begin(directory, root):
            directory.mkdir(parents=True)
            return {"git_commit": "a" * 40}
        def subprocess_run(command, **kwargs):
            if "scripts/taxonomy/run_taxonomy.py" in command:
                target = command[command.index("--output-dir") + 1]
                (self.root / target).mkdir(parents=True)
                events.append("classify_" + command[command.index("--threads") + 1])
            elif "scripts/taxonomy/audit_repetition.py" in command:
                self.assertIn("--require-identical", command)
                events.append("compare")
        def native(context, root, repeat=False, expected_revision=None):
            self.assertEqual(expected_revision, "a" * 40)
            events.append("validate_repeat" if repeat else "validate_primary")
            return self.native
        def finish(*args):
            events.append("finish")
            (audit / "SUCCESS").write_text("a" * 40)
        with patch.object(TAX, "check_sources", return_value={"status": "valid"}), patch.object(TAX, "inputs", return_value=[self.policy_path]), \
             patch.object(TAX, "begin_stage", side_effect=begin), patch.object(TAX.subprocess, "run", side_effect=subprocess_run), \
             patch.object(TAX, "native_result", side_effect=native), patch.object(TAX, "check_repetition", return_value={"all_identical": True}), \
             patch.object(TAX, "collect_versions", return_value={"Python": "fixture"}), patch.object(TAX, "finish_stage", side_effect=finish), \
             patch.object(TAX, "validate", side_effect=ValueError("Final validation failed")):
            with self.assertRaisesRegex(ValueError, "Final validation failed"):
                TAX.execute(self.policy_path, self.root)
        self.assertEqual(events, ["classify_2", "validate_primary", "classify_1", "validate_repeat", "compare", "finish"])
        self.assertFalse((audit / "SUCCESS").exists())
        self.assertFalse((audit / "validation.json").exists())
        self.assertTrue((audit / "failure.json").exists())


if __name__ == "__main__":
    unittest.main()
