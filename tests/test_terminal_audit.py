"""Synthetic terminal-trimming fixtures are structural tests, not study evidence."""

import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("terminal_audit", ROOT / "scripts/terminal_qc/audit.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def fastq(header, sequence, quality=None):
    if quality is None:
        quality = "I" * len(sequence)
    return f"@{header}\n{sequence}\n+\n{quality}\n"


def fastqc_text():
    columns = "\t".join(AUDIT.QUALITY_COLUMNS)
    return (
        "##FastQC\t0.12.1\n"
        ">>Basic Statistics\tpass\n#Measure\tValue\n"
        "Filename\tsynthetic.fastq.gz\nTotal Sequences\t4\nSequence length\t0-5\n>>END_MODULE\n"
        f">>Per base sequence quality\tfail\n#{columns}\n"
        "1\t30\t31\t28\t34\t20\t36\n"
        "2-4\t20\t21\t18\t24\t10\t26\n"
        ">>END_MODULE\n>>Adapter Content\twarn\n#Position\tSynthetic adapter\n"
        "1\t0\n>>END_MODULE\n"
    )


class PairedTerminalAuditTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        root = Path(self.directory.name)
        self.inputs = [root / "before_1.fastq.gz", root / "before_2.fastq.gz"]
        self.outputs = [root / "after_1.fastq.gz", root / "after_2.fastq.gz"]
        self.records = [
            fastq("a/1 extra", "AACCGG") + fastq("b/1", "GG") + fastq("c/1", "AAAA"),
            fastq("a/2 extra", "TTTCCC") + fastq("b/2", "CCC") + fastq("c/2", "TTTT"),
            fastq("a/1 extra", "AACC") + fastq("b/1", "") + fastq("c/1", "AAAA"),
            fastq("a/2 extra", "TTT") + fastq("b/2", "") + fastq("c/2", "TTTT"),
        ]
        self.write_records()

    def write_records(self):
        for path, records in zip([*self.inputs, *self.outputs], self.records):
            with gzip.open(path, "wt", encoding="ascii", newline="\n") as handle:
                handle.write(records)

    def audit(self):
        return AUDIT.audit_pair(self.inputs, self.outputs, [2, 3], min_len=4)

    def test_exact_terminal_removal_empty_reads_and_threshold_crossing(self):
        result = self.audit()
        self.assertEqual(result["pairs"], 3)
        self.assertEqual(result["pairs_below_min_len_baseline"], 1)
        self.assertEqual(result["pairs_below_min_len"], 2)
        self.assertEqual(result["pairs_newly_below_min_len"], 1)
        r1, r2 = result["directions"]["R1"], result["directions"]["R2"]
        self.assertEqual(r1["changed_reads"], 2)
        self.assertEqual(r1["bases_removed"], 4)
        self.assertEqual(r1["length_counts"], {"0": 1, "4": 2})
        self.assertEqual(r1["baseline_length_counts"], {"2": 1, "4": 1, "6": 1})
        self.assertEqual((r1["min_length"], r1["max_length"]), (0, 4))
        self.assertEqual((r1["baseline_min_length"], r1["baseline_max_length"]), (2, 6))
        self.assertEqual(r1["baseline_below_min_len"], 1)
        self.assertEqual(r1["newly_below_min_len"], 0)
        self.assertEqual(r2["below_min_len"], 2)
        self.assertEqual(r2["newly_below_min_len"], 1)
        self.assertEqual(r2["bases_removed"], 6)

    def test_casava_whitespace_mate_identifiers_are_normalized(self):
        self.records = [
            fastq("cluster 1:N:0:TAG", "AA"), fastq("cluster 2:N:0:TAG", "TTT"),
            fastq("cluster 1:N:0:TAG", "AA"), fastq("cluster 2:N:0:TAG", "TTT"),
        ]
        self.write_records()
        self.assertEqual(self.audit()["pairs"], 1)

    def test_mutated_retained_sequence_or_quality_fails(self):
        for replacement in (fastq("a/1 extra", "AGCC"), fastq("a/1 extra", "AACC", "IIHI")):
            with self.subTest(replacement=replacement):
                self.records[2] = replacement + fastq("b/1", "") + fastq("c/1", "AAAA")
                self.write_records()
                with self.assertRaisesRegex(ValueError, "altered retained"):
                    self.audit()

    def test_changed_header_or_separator_fails(self):
        for old, new in (("a/1 extra", "a/1 changed"), ("\n+\n", "\n+changed\n")):
            with self.subTest(change=new):
                original = self.records[2]
                self.records[2] = original.replace(old, new, 1)
                self.write_records()
                with self.assertRaisesRegex(ValueError, "header or separator changed"):
                    self.audit()
                self.records[2] = original

    def test_non_primer_length_removal_fails(self):
        self.records[2] = self.records[2].replace("AACC\n+\nIIII", "AAC\n+\nIII")
        self.write_records()
        with self.assertRaisesRegex(ValueError, "Unexpected terminal trimming length"):
            self.audit()

    def test_desynchronized_input_mates_fail(self):
        self.records[1] = self.records[1].replace("a/2", "other/2", 1)
        self.write_records()
        with self.assertRaisesRegex(ValueError, "different identifiers or order"):
            self.audit()

    def test_swapped_mate_suffix_fails(self):
        self.records[1] = self.records[1].replace("a/2", "a/1", 1)
        self.write_records()
        with self.assertRaisesRegex(ValueError, "mate suffix"):
            self.audit()

    def test_lost_read_fails(self):
        self.records[3] = self.records[3].removesuffix(fastq("c/2", "TTTT"))
        self.write_records()
        with self.assertRaisesRegex(ValueError, "different numbers"):
            self.audit()

    def test_bad_fastq_structure_fails(self):
        invalid_records = (
            "a\nA\n+\nI\n", "@a\nA\n-\nI\n", "@a\nAA\n+\nI\n",
            "@a\nA\n+\n", "@a\nA\n+\n \n", "@\nA\n+\nI\n",
        )
        for value in invalid_records:
            with self.subTest(value=value):
                self.records[0] = value
                self.write_records()
                with self.assertRaises(ValueError):
                    self.audit()

    def test_truncated_gzip_fails_instead_of_returning_partial_counts(self):
        path = self.inputs[0]
        path.write_bytes(path.read_bytes()[:-5])
        with self.assertRaisesRegex(ValueError, "complete gzip"):
            self.audit()

    def test_empty_streams_have_no_invented_minimum_length(self):
        self.records = ["", "", "", ""]
        self.write_records()
        result = self.audit()
        self.assertEqual(result["pairs"], 0)
        self.assertIsNone(result["directions"]["R1"]["min_length"])
        self.assertEqual(result["directions"]["R1"]["length_counts"], {})


class FastqcCoverageTests(unittest.TestCase):
    def parse(self, content, duplicate=False):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.zip"
            with ZipFile(path, "w") as archive:
                archive.writestr("synthetic/fastqc_data.txt", content)
                if duplicate:
                    archive.writestr("another/fastqc_data.txt", content)
            return AUDIT.read_fastqc(path)

    def test_complete_report_preserves_native_columns_and_failures(self):
        result = self.parse(fastqc_text())
        self.assertEqual(result["basic"]["Total Sequences"], "4")
        self.assertEqual(result["module_status"]["Per base sequence quality"], "fail")
        self.assertEqual(result["module_status"]["Adapter Content"], "warn")
        self.assertEqual(result["quality_rows"][1]["Base"], "2-4")
        self.assertEqual(result["quality_rows"][1]["10th Percentile"], "10")

    def test_multiple_data_members_fail(self):
        with self.assertRaisesRegex(ValueError, "exactly one"):
            self.parse(fastqc_text(), duplicate=True)

    def test_incomplete_modules_wrong_headers_and_total_reads_fail(self):
        original = fastqc_text()
        for invalid in (
            original.removesuffix(">>END_MODULE\n"),
            original.replace("#Measure\tValue", "#Measure\tDifferent"),
            original.replace("#Base\tMean", "#Position\tMean"),
            original.replace("Total Sequences\t4", "Total Sequences\t-4"),
            original.replace("Total Sequences\t4", "Total Sequences\t4.5"),
            original.replace("##FastQC\t0.12.1\n", ""),
            original.replace("1\t30\t31", "1\tInfinity\t31"),
        ):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                self.parse(invalid)

    def test_native_missing_percentiles_are_retained_without_imputation(self):
        content = fastqc_text().replace("1\t30\t31\t28\t34\t20\t36", "1\t30\tNaN\tNaN\tNaN\tNaN\tNaN")
        row = self.parse(content)["quality_rows"][0]
        self.assertEqual(row["Mean"], "30")
        self.assertEqual(row["Median"], "NaN")
        self.assertEqual(row["90th Percentile"], "NaN")

    def test_variable_length_coverage_counts_each_observed_base(self):
        rows = [{"Base": "1", "Mean": "31"}, {"Base": "2-4", "Mean": "23"}, {"Base": "5-6", "Mean": "10"}]
        observed = AUDIT.quality_coverage(rows, {"0": 1, "1": 1, "3": 1, "5": 1})
        self.assertEqual(observed[0]["reads_reaching_start"], 3)
        self.assertEqual(observed[0]["fraction_reads_reaching_start"], 0.75)
        self.assertEqual(observed[1]["reads_reaching_start"], 2)
        self.assertEqual(observed[1]["reads_reaching_end"], 1)
        self.assertEqual(observed[1]["base_observations"], 5)
        self.assertEqual(observed[1]["fraction_reads_reaching_start"], 0.5)
        self.assertEqual(observed[2]["reads_reaching_end"], 0)
        self.assertEqual(observed[2]["base_observations"], 1)
        self.assertEqual(observed[1]["Mean"], "23")
        self.assertEqual(rows[1], {"Base": "2-4", "Mean": "23"})

    def test_no_observed_bases_requires_an_empty_quality_report(self):
        for histogram in ({}, {"0": 4}, {"5": 0}, {"0": 0, "5": 0}):
            with self.subTest(histogram=histogram):
                self.assertEqual(AUDIT.quality_coverage([], histogram), [])
                with self.assertRaisesRegex(ValueError, "no observed bases"):
                    AUDIT.quality_coverage([{"Base": "1-5"}], histogram)

    def test_missing_initial_or_internal_quality_bins_fail(self):
        for rows in ([{"Base": "2-5"}], [{"Base": "1"}, {"Base": "3-5"}]):
            with self.subTest(rows=rows), self.assertRaisesRegex(ValueError, "contiguous from base 1"):
                AUDIT.quality_coverage(rows, {"5": 4})
        for text in (
            fastqc_text().replace("1\t30\t31", "2\t30\t31"),
            fastqc_text().replace("2-4\t20", "3-4\t20"),
        ):
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, "contiguous from base 1"):
                self.parse(text)

    def test_empty_or_truncated_quality_reports_cannot_omit_observed_bases(self):
        for rows in ([], [{"Base": "1"}], [{"Base": "1"}, {"Base": "2-4"}]):
            with self.subTest(rows=rows), self.assertRaisesRegex(ValueError, "omits observed terminal bases"):
                AUDIT.quality_coverage(rows, {"0": 1, "5": 3})
        with self.assertRaisesRegex(ValueError, "no observed bases"):
            AUDIT.quality_coverage([{"Base": "1-5"}, {"Base": "6-10"}], {"5": 4})

    def test_invalid_bins_or_histograms_fail(self):
        for rows, counts in (
            ([{"Base": "0"}], {"1": 1}), ([{"Base": "4-2"}], {"1": 1}),
            ([{"Base": "1-3"}, {"Base": "3-4"}], {"1": 1}),
            ([{"Base": "1"}], {"-1": 1}), ([{"Base": "1"}], {"1": 1.5}),
            ([{"Base": "1"}], {"1": -1}), ([{"Base": "1"}], {"1": 1, 1: 2}),
        ):
            with self.subTest(rows=rows, counts=counts), self.assertRaises(ValueError):
                AUDIT.quality_coverage(rows, counts)


if __name__ == "__main__":
    unittest.main()
