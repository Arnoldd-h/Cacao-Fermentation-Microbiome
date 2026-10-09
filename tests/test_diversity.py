"""Synthetic fixtures check conservation and analytical formulas, never real evidence."""

import copy
import json
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts/filtering"))
from validate_bacterial_table import validate_tables
from cacao_inventory.analysis_artifacts import finish_stage, read_counts, read_table, validate_artifacts, validate_policy, write_json
from cacao_inventory.bacterial_filter import prepare_tables
from cacao_inventory.config import load_json_yaml
from cacao_inventory.diversity_validation import validate_diversity_tables
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.provenance import file_records


def fixture(root):
    config = load_json_yaml(ROOT / "config/diversity.yaml")
    paths = validate_policy(config, root)
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    ranks = config["filtering"]["tax_levels"]
    features = ["a", "b", "c", "org", "arch", "unknown"]
    values = [[5, 0, 1, 3, 1, 0], [0, 4, 1, 1, 0, 2], [2, 1, 2, 0, 0, 1]]
    rows = [{"study_id": "synthetic", "sample_id": f"s{i}", **dict(zip(features, map(str, vector)))} for i, vector in enumerate(values)]
    seq = [{"study_id": "synthetic", "asv_id": f, "sequence": "ACGT" * (i + 1), "total_reads": str(sum(v[i] for v in values))} for i, f in enumerate(features)]
    taxa = []
    for feature in features:
        lineage = dict(zip(ranks, ["Bacteria", "P", "C", "O", "F", "G_" + feature]))
        if feature == "org":
            lineage.update(Order="Chloroplast", Family="", Genus="")
        if feature == "arch":
            lineage["Kingdom"] = "Archaea"
        if feature == "unknown":
            lineage = dict.fromkeys(ranks, "")
        taxa.append({"study_id": "synthetic", "asv_id": feature, **lineage})
    metadata = [{"study_id": "synthetic", "sample_id": f"s{i}", "run_accession": f"fixture_run_{i}",
                 "fermentation_batch": f"fixture_batch_{i}", "fermentation_hours": str(i * 24),
                 "relative_time": str(i / 2), "fermentation_stage": stage, "sampling_stratum": ""}
                for i, stage in enumerate(("early", "mid", "late"))]
    for directory, name, table in ((paths["dada2_dir"], "asv_counts.tsv", rows),
                                   (paths["dada2_dir"], "asv_sequences.tsv", seq),
                                   (paths["dada2_dir"], "sample_metadata.tsv", metadata),
                                   (paths["taxonomy_dir"], "taxonomy.tsv", taxa)):
        write_tsv_atomic(directory / name, table, list(table[0]))
    return config, paths


def materialize(config, paths):
    tables, summary = prepare_tables(paths["dada2_dir"], paths["taxonomy_dir"], config["filtering"])
    for name, rows in tables.items():
        target = paths["filtering_dir"] if name in ("asv_filter_log.tsv", "sample_retention.tsv") else paths["bacterial_dir"]
        write_tsv_atomic(target / name, rows, list(rows[0]))
    write_json(paths["filtering_dir"] / "summary.json", summary)
    return tables, summary


class BacterialSeparationTests(unittest.TestCase):
    def test_organelles_and_known_non_bacteria_are_removed_unknowns_and_rare_asvs_retained(self):
        with tempfile.TemporaryDirectory() as temporary:
            config, paths = fixture(Path(temporary))
            tables, summary = materialize(config, paths)
            self.assertEqual(summary["excluded_asvs"], 2)
            self.assertEqual(summary["input_reads"], summary["retained_reads"] + summary["excluded_reads"])
            self.assertEqual([row["asv_id"] for row in tables["asv_sequences.tsv"]], ["a", "b", "c", "unknown"])
            self.assertEqual(summary["unclassified_kingdom_retained_asvs"], 1)
            self.assertEqual(validate_tables(config, paths), summary)

    def test_empty_sample_fails_without_silently_removing_it(self):
        with tempfile.TemporaryDirectory() as temporary:
            config, paths = fixture(Path(temporary))
            path = paths["dada2_dir"] / "asv_counts.tsv"
            rows, features = read_counts(path)
            for f in ("a", "c"):
                rows[1][f] = str(int(rows[1][f]) + int(rows[0][f]))
                rows[0][f] = "0"
            write_tsv_atomic(path, rows, list(rows[0]))
            with self.assertRaisesRegex(ValueError, "empties sample"):
                prepare_tables(paths["dada2_dir"], paths["taxonomy_dir"], config["filtering"])

    def test_count_or_taxonomy_tampering_fails_independent_validation(self):
        for name in ("asv_counts.tsv", "taxonomy.tsv"):
            with tempfile.TemporaryDirectory() as temporary:
                config, paths = fixture(Path(temporary))
                materialize(config, paths)
                path = paths["bacterial_dir"] / name
                rows = read_table(path, [])
                rows[0]["a" if name == "asv_counts.tsv" else "Genus"] = "99" if name == "asv_counts.tsv" else "invented"
                write_tsv_atomic(path, rows, list(rows[0]))
                with self.assertRaises(ValueError):
                    validate_tables(config, paths)

    def test_sample_metadata_cannot_be_invented(self):
        with tempfile.TemporaryDirectory() as temporary:
            config, paths = fixture(Path(temporary))
            materialize(config, paths)
            path = paths["bacterial_dir"] / "sample_metadata.tsv"
            rows = read_table(path, [])
            rows[0]["fermentation_batch"] = "invented"
            write_tsv_atomic(path, rows, list(rows[0]))
            with self.assertRaisesRegex(ValueError, "metadata changed"):
                validate_tables(config, paths)

    def test_partial_name_match_does_not_exclude_a_taxon(self):
        with tempfile.TemporaryDirectory() as temporary:
            config, paths = fixture(Path(temporary))
            path = paths["taxonomy_dir"] / "taxonomy.tsv"
            rows = read_table(path, [])
            rows[0]["Genus"] = "Chloroplast_like_fixture"
            write_tsv_atomic(path, rows, list(rows[0]))
            tables, _ = materialize(config, paths)
            self.assertIn("a", tables["asv_counts.tsv"][0])

    def test_mixed_studies_and_invalid_count_headers_fail(self):
        for text in ("study_id\tsample_id\ta\ta\nx\ts\t1\t2\n",
                     "study_id\tsample_id\ta\nx\ts\t1.5\n",
                     "study_id\tsample_id\ta\nx\ts\t1\t2\n",
                     "study_id\tsample_id\ta\nx\ts\t1\ny\tt\t2\n"):
            with tempfile.TemporaryDirectory() as temporary:
                path = Path(temporary) / "counts.tsv"
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_counts(path)

    def test_inference_zero_constant_and_raw_output_path_are_rejected(self):
        for section, field, value in (("analysis", "inference", True), ("analysis", "primary_pseudocount", 0),
                                      ("paths", "bacterial_dir", "data/raw/fixture"),
                                      ("analysis", "sensitivity_pseudocounts", [1])):
            with tempfile.TemporaryDirectory() as temporary:
                config, _ = fixture(Path(temporary))
                config[section][field] = value
                with self.assertRaises(ValueError):
                    validate_policy(config, Path(temporary))


class AnalysisProvenanceTests(unittest.TestCase):
    def test_checksum_change_and_success_mismatch_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = root / "results/fixture"
            report.mkdir(parents=True)
            source = root / "source.txt"
            source.write_text("synthetic input\n")
            config_path = root / "config.json"
            config = {"synthetic_fixture": True}
            write_json(config_path, config)
            output = report / "data.tsv"
            output.write_text("synthetic output\n")
            before = file_records([source, config_path], root)
            state = {"git_commit": "a" * 40, "git_dirty": False, "git_status_at_start": [], "started_at_utc": "fixture"}
            finish_stage(report, root, state, before, config_path, config, [output], {}, {}, [])
            expected = [output, report / "config_snapshot.yaml", report / "input_checksums.json"]
            self.assertEqual(validate_artifacts(report, root, expected, config)["status"], "valid")
            (report / "SUCCESS").write_text("b" * 40)
            with self.assertRaisesRegex(ValueError, "success"):
                validate_artifacts(report, root, expected, config)
            (report / "SUCCESS").write_text("a" * 40)
            output.write_text("tampered\n")
            with self.assertRaisesRegex(ValueError, "checksum"):
                validate_artifacts(report, root, expected, config)

    def test_changed_input_during_execution_cannot_create_success(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "config.json"
            write_json(source, {"synthetic": True})
            before = file_records([source], root)
            source.write_text("changed\n")
            with self.assertRaisesRegex(ValueError, "Inputs changed"):
                finish_stage(root / "report", root, {}, before, source, {}, [], {}, {}, [])
            self.assertFalse((root / "report/SUCCESS").exists())


@unittest.skipUnless(shutil.which("Rscript"), "R integration requires the declared Linux environment")
class DiversityFormulaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="diversity fixture root with spaces ")
        cls.root = Path(cls.temporary.name)
        cls.config, cls.paths = fixture(cls.root)
        materialize(cls.config, cls.paths)
        # Three disjoint, equally abundant synthetic communities have known
        # entropy ln(2), Gini-Simpson 1/2 and Hill numbers 2.
        rows = [{"study_id": "synthetic", "sample_id": f"s{i}", "a": str(v[0]), "b": str(v[1]), "c": str(v[2])}
                for i, v in enumerate(((1, 1, 0), (0, 1, 1), (1, 0, 1)))]
        write_tsv_atomic(cls.paths["bacterial_dir"] / "asv_counts.tsv", rows, list(rows[0]))
        config_path = cls.root / "config/diversity.yaml"
        write_json(config_path, cls.config)
        for name in ("scripts/diversity/describe_pilot.R", "scripts/dada2/helpers.R"):
            destination = cls.root / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / name, destination)
        result = subprocess.run(["Rscript", str(cls.root / "scripts/diversity/describe_pilot.R"), "--config", str(config_path)],
                                capture_output=True, text=True, timeout=90)
        if result.returncode:
            cls.temporary.cleanup()
            raise AssertionError(result.stdout + result.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def validate(self):
        return validate_diversity_tables(self.paths["diversity_dir"], self.paths["bacterial_dir"], self.config)

    def test_r_metrics_and_geometry_match_known_communities(self):
        summary = self.validate()
        self.assertEqual(summary["retained_reads"], 6)
        alpha = read_table(self.paths["diversity_dir"] / "alpha_diversity.tsv", [])
        self.assertAlmostEqual(float(alpha[0]["shannon"]), math.log(2))
        self.assertAlmostEqual(float(alpha[0]["inverse_simpson"]), 2)
        distance = read_table(self.paths["diversity_dir"] / "aitchison_distances.tsv", [])
        pair = next(row for row in distance if row["sample_id_1"] == "s0" and row["sample_id_2"] == "s1" and float(row["pseudocount"]) == 1)
        self.assertAlmostEqual(float(pair["distance"]), math.sqrt(2) * math.log(2))

    def test_changed_alpha_clr_distances_and_pca_fail_without_trusting_r_summaries(self):
        for name, field in (("alpha_diversity.tsv", "shannon"), ("clr_coordinates.tsv", "clr"),
                            ("aitchison_distances.tsv", "distance"), ("bray_curtis_distances.tsv", "distance"),
                            ("pca_scores.tsv", "score"), ("pca_variance.tsv", "explained_fraction")):
            path = self.paths["diversity_dir"] / name
            original = path.read_bytes()
            try:
                rows = read_table(path, [])
                rows[0][field] = str(float(rows[0][field]) + 1)
                write_tsv_atomic(path, rows, list(rows[0]))
                with self.assertRaises(ValueError):
                    self.validate()
            finally:
                path.write_bytes(original)

    def test_missing_pseudocount_or_duplicate_sample_fails(self):
        for name in ("clr_coordinates.tsv", "alpha_diversity.tsv"):
            path = self.paths["diversity_dir"] / name
            original = path.read_bytes()
            try:
                rows = read_table(path, [])
                rows = [row for row in rows if float(row["pseudocount"]) == 1] if name.startswith("clr") else [*rows, rows[0]]
                write_tsv_atomic(path, rows, list(rows[0]))
                with self.assertRaises(ValueError):
                    self.validate()
            finally:
                path.write_bytes(original)


if __name__ == "__main__":
    unittest.main()
