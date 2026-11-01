"""Reject full-study QC corruption using copies of validated pilot summaries."""

import copy
import gzip
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.analysis_artifacts import read_table
from cacao_inventory.config import load_json_yaml
from cacao_inventory.io import write_tsv_atomic
from cacao_inventory.pilot_workflow import load_pilot_manifest
from cacao_inventory.qc import cutadapt_command
from cacao_inventory.study_qc import cutadapt_rows, reconcile_qc, trimmed_path


class StudyQcTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.config = load_json_yaml(ROOT / "config/full_study.yaml")
        self.project = load_json_yaml(ROOT / self.config["project_config"])
        self.manifest = load_pilot_manifest(ROOT / "metadata/pilot_manifest.tsv")[:2]
        self.runs = {row["run_accession"] for row in self.manifest}
        self.directory = self.root / self.config["qc_dir"]
        self.tables = {}
        for name in ("raw_read_quality", "trimmed_read_quality", "cutadapt_summary", "primer_detection_raw", "primer_detection_trimmed"):
            source = "primer_detection" if name == "primer_detection_raw" else name
            self.tables[name] = [row for row in read_table(ROOT / "results/qc/pilot" / (source + ".tsv"), []) if row["run_accession"] in self.runs]
            self.write(name)
        parser = patch("cacao_inventory.study_qc.cutadapt_rows", return_value=copy.deepcopy(self.tables["cutadapt_summary"]))
        parser.start()
        self.addCleanup(parser.stop)

    def write(self, name):
        write_tsv_atomic(self.directory / (name + ".tsv"), self.tables[name], list(self.tables[name][0]))

    def test_consistent_qc_conserves_all_eligible_runs_and_pairs(self):
        summary = reconcile_qc(self.config, self.project, self.manifest, self.root)
        expected = sum(int(row["total_sequences"]) for row in self.tables["raw_read_quality"] if row["read_direction"] == "R1")
        self.assertEqual(summary["input_pairs"], expected)
        self.assertEqual(summary["trimmed_pairs"], expected)
        self.assertEqual(summary["runs"], 2)

    def test_altered_counts_or_sample_metadata_fail(self):
        for field, value in (("total_sequences", "0"), ("sample_id", "invented")):
            original = self.tables["trimmed_read_quality"][0][field]
            self.tables["trimmed_read_quality"][0][field] = value
            self.write("trimmed_read_quality")
            with self.subTest(field=field), self.assertRaises(ValueError):
                reconcile_qc(self.config, self.project, self.manifest, self.root)
            self.tables["trimmed_read_quality"][0][field] = original

    def test_missing_run_fails(self):
        self.tables["raw_read_quality"].pop()
        self.write("raw_read_quality")
        with self.assertRaises(ValueError):
            reconcile_qc(self.config, self.project, self.manifest, self.root)

    def test_altered_cutadapt_parameters_fail_independent_reconciliation(self):
        self.tables["cutadapt_summary"][0]["error_rate"] = "0.5"
        self.write("cutadapt_summary")
        with self.assertRaises(ValueError):
            reconcile_qc(self.config, self.project, self.manifest, self.root)

    def test_altered_screen_design_counts_and_parameters_fail(self):
        for field, value in (("primer", "unknown"), ("reads_examined", "1"), ("error_rate", "0.5"), ("match_percent", "-2"), ("sample_id", "invented")):
            original = self.tables["primer_detection_trimmed"][0][field]
            self.tables["primer_detection_trimmed"][0][field] = value
            self.write("primer_detection_trimmed")
            with self.subTest(field=field), self.assertRaises(ValueError):
                reconcile_qc(self.config, self.project, self.manifest, self.root)
            self.tables["primer_detection_trimmed"][0][field] = original
        self.tables["primer_detection_trimmed"][0] = self.tables["primer_detection_trimmed"][1].copy()
        self.write("primer_detection_trimmed")
        with self.assertRaises(ValueError):
            reconcile_qc(self.config, self.project, self.manifest, self.root)


@unittest.skipUnless(shutil.which("cutadapt"), "Requires the declared Cutadapt environment")
class StudyCutadaptRoutingTests(unittest.TestCase):
    def test_real_cutadapt_report_validates_in_registered_full_study_directory(self):
        """Tiny synthetic reads exercise the unmocked full-study parser wiring."""
        config = load_json_yaml(ROOT / "config/full_study.yaml")
        project = load_json_yaml(ROOT / config["project_config"])
        row = load_pilot_manifest(ROOT / "metadata/pilot_manifest.tsv")[0]
        processing = project["amplicon_processing"][row["bioproject"]]
        choices = {"R": "A", "Y": "C", "S": "C", "W": "A", "K": "G", "M": "A", "B": "C", "D": "A", "H": "A", "V": "A", "N": "A"}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            outputs = {}
            for direction, suffix, role in (("R1", "1", "forward"), ("R2", "2", "reverse")):
                primer = processing["primers"][role]["gene_specific_sequence"]
                sequence = "".join(choices.get(base, base) for base in primer) + "ACGT" * 30
                other = "T" * 120
                reads = f"@synthetic_1/{suffix}\n{sequence}\n+\n{'I' * len(sequence)}\n@synthetic_2/{suffix}\n{other}\n+\n{'I' * len(other)}\n"
                source = root / "data/raw" / row["study_id"] / row["run_accession"] / (row["run_accession"] + "_" + suffix + ".fastq.gz")
                source.parent.mkdir(parents=True, exist_ok=True)
                source.write_bytes(gzip.compress(reads.encode(), mtime=0))
                output = Path(trimmed_path(config, row, suffix))
                (root / output).parent.mkdir(parents=True, exist_ok=True)
                outputs[direction] = output.as_posix()
            report = Path(config["qc_dir"]) / "cutadapt" / (row["run_accession"] + ".cutadapt.json")
            (root / report).parent.mkdir(parents=True)
            command = cutadapt_command(row, processing["primers"], processing["cutadapt"], outputs, report.as_posix())
            result = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            summary = cutadapt_rows(config, project, [row], root)
            self.assertEqual(len(summary), 2)
            self.assertEqual({item["output_file"] for item in summary}, set(outputs.values()))
            self.assertTrue(all(item["input_reads"] == item["reads_written"] == "2" for item in summary))
            self.assertTrue(all(item["reads_with_adapter"] == "1" for item in summary))


if __name__ == "__main__":
    unittest.main()
