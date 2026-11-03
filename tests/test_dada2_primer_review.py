"""Synthetic primer motifs and independent native/outer audit fixtures."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("dada2_primer_review", ROOT / "scripts/dada2_primer_review/scan.py")
SCAN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SCAN)


def processing():
    return {"study_id": "synthetic", "primers": {
        "forward": {"name": "forward", "gene_specific_sequence": "ACGN"},
        "reverse": {"name": "reverse", "gene_specific_sequence": "ACGT"}},
        "primer_detection": {"error_rate": 0, "allow_indels": False, "minimum_overlap_fraction": 1}}


def tables():
    sequences = [{"study_id": "synthetic", "asv_id": "A", "sequence": "ACGTACGT", "length": "8", "total_reads": "5"},
                 {"study_id": "synthetic", "asv_id": "B", "sequence": "GGGG", "length": "4", "total_reads": "5"}]
    counts = [{"study_id": "synthetic", "sample_id": "s1", "A": "2", "B": "5"},
              {"study_id": "synthetic", "sample_id": "s2", "A": "3", "B": "0"}]
    return sequences, counts


class PrimerScannerTests(unittest.TestCase):
    def test_all_overlapping_full_positions_without_partial_matches(self):
        self.assertEqual(SCAN.find_matches("AAAAAA", "AAAA", 0),
                         [{"start": n, "end": n + 4, "errors": 0} for n in (0, 1, 2)])
        self.assertEqual(SCAN.find_matches("AAA", "AAAA", 0), [])
        self.assertEqual(SCAN.find_matches("AAAACCCC", "AACCC", 0), [{"start": 2, "end": 7, "errors": 0}])

    def test_adapter_iupac_and_n_effective_error_denominator(self):
        self.assertEqual(SCAN.find_matches("AAGCT", "RRSSW", 0), [{"start": 0, "end": 5, "errors": 0}])
        self.assertEqual(SCAN.maximum_errors("ACGTNNNNNN", 0.1), 0)
        self.assertEqual(SCAN.maximum_errors("ACGTRYSWKM", 0.1), 1)
        self.assertEqual(SCAN.maximum_errors("CCTACGGGNGGCWGCAG", 0.1), 1)
        self.assertEqual(SCAN.maximum_errors("GACTACHVGGGTATCTAATCC", 0.1), 2)
        self.assertEqual(SCAN.find_matches("TTTT", "ANNN", 0.1), [])

    def test_invalid_sequences_patterns_or_rates_fail(self):
        for sequence, pattern, rate in (("ACNT", "ACGT", 0), ("acgt", "ACGT", 0), ("", "ACGT", 0),
                                        ("ACGT", "NNNN", 0.1), ("ACGT", "ACGX", 0),
                                        ("ACGT", "ACGT", True), ("ACGT", "ACGT", "nan"),
                                        ("ACGT", "ACGT", -0.1), ("ACGT", "ACGT", "bad")):
            with self.subTest(sequence=sequence, pattern=pattern, rate=rate), self.assertRaises(ValueError):
                SCAN.find_matches(sequence, pattern, rate)

    def test_read_and_asv_unions_do_not_double_count_repeated_or_overlapping_patterns(self):
        coordinates, samples, summary = SCAN.scan_tables(*tables(), processing())
        self.assertEqual(len(coordinates), 8)
        self.assertEqual(summary["matched_asvs"], 1)
        self.assertEqual(summary["matched_reads"], 5)
        self.assertEqual(summary["total_reads"], 10)
        self.assertEqual([row["matched_reads"] for row in samples], [2, 3])
        self.assertEqual([row["matched_asvs"] for row in samples], [1, 1])
        self.assertEqual([row["matched_reads"] for row in summary["patterns"]], [5] * 4)
        self.assertEqual({row["orientation"] for row in coordinates}, {"as_synthesized", "reverse_complement"})

    def test_no_matches_keeps_every_sample_and_count(self):
        sequences, counts = tables()
        sequences[0]["sequence"] = "CCCCCCCC"
        coordinates, samples, summary = SCAN.scan_tables(sequences, counts, processing())
        self.assertEqual(coordinates, [])
        self.assertEqual(len(samples), 2)
        self.assertEqual(summary["matched_reads"], 0)
        self.assertEqual(summary["total_reads"], 10)
        self.assertEqual(summary["asv_exclusions"], 0)

    def test_inconsistent_counts_ids_sequences_or_profile_fail(self):
        for kind in ("duplicate_asv", "duplicate_sample", "missing_asv", "changed_total", "ambiguous_asv", "mixed_study", "indels", "partial"):
            sequences, counts = tables()
            profile = processing()
            if kind == "duplicate_asv": sequences.append(deepcopy(sequences[0]))
            elif kind == "duplicate_sample": counts.append(deepcopy(counts[0]))
            elif kind == "missing_asv": del counts[0]["A"]
            elif kind == "changed_total": counts[0]["A"] = "3"
            elif kind == "ambiguous_asv": sequences[0]["sequence"] = "ACGNACGT"
            elif kind == "mixed_study": counts[0]["study_id"] = "other"
            elif kind == "indels": profile["primer_detection"]["allow_indels"] = True
            else: profile["primer_detection"]["minimum_overlap_fraction"] = 0.5
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                SCAN.scan_tables(sequences, counts, profile)

    @unittest.skipUnless(shutil.which("cutadapt"), "Cutadapt integration requires the registered environment")
    def test_hamming_thresholds_and_iupac_positions_match_cutadapt(self):
        gene_primers = ("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC", "ACGTNNNNNN")
        probes = [pattern for primer in gene_primers for pattern in (primer, primer.translate(SCAN.COMPLEMENT)[::-1])]
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for pattern in probes:
                concrete = "".join(SCAN.IUPAC[base][0] for base in pattern)
                informative = [i for i, base in enumerate(pattern) if base != "N"]
                sequences = []
                for errors in range(SCAN.maximum_errors(pattern, 0.1) + 2):
                    bases = list(concrete)
                    for offset in informative[:errors]:
                        bases[offset] = next(base for base in "ACGT" if base not in SCAN.IUPAC[pattern[offset]])
                    sequences.append("TTTTTC" + "".join(bases) + "GGGGGG")
                path, info = directory / "synthetic.fastq", directory / "matches.tsv"
                path.write_text("".join(f"@case_{i}\n{sequence}\n+\n{'I' * len(sequence)}\n" for i, sequence in enumerate(sequences)))
                subprocess.run(["cutadapt", "--action=none", "--no-indels", "--error-rate", "0.1", "--overlap", str(len(pattern)),
                                "--front", "probe=" + pattern, "--info-file", str(info), "-o", str(directory / "unchanged.fastq"), str(path)],
                               check=True, capture_output=True, text=True)
                for index, line in enumerate(info.read_text().splitlines()):
                    fields = line.split("\t")
                    found = SCAN.find_matches(sequences[index], pattern, 0.1)
                    observed = [] if fields[1] == "-1" else [{"start": int(fields[2]), "end": int(fields[3]), "errors": int(fields[1])}]
                    with self.subTest(pattern=pattern, sequence=sequences[index]):
                        self.assertEqual(found, observed)


class SourceAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.directory = self.root / "results/dada2/synthetic/run"
        self.audit = self.root / "results/dada2/synthetic/audit"
        self.directory.mkdir(parents=True)
        self.audit.mkdir()
        self.revision = "a" * 40
        self.policy = {"output_dir": "results/dada2/synthetic/run", "audit_dir": "results/dada2/synthetic/audit"}
        self.scope = {"study_id": "synthetic", "expected_runs": 1}
        self.source = self.root / "config/source.json"
        self.source.parent.mkdir()
        self.source.write_text("{}\n")
        native_path = self.root / "scripts/dada2/validate_outputs.py"
        native_path.parent.mkdir(parents=True)
        shutil.copyfile(ROOT / "scripts/dada2/validate_outputs.py", native_path)
        content = {
            "asv_counts.tsv": "study_id\tsample_id\tA\nsynthetic\ts1\t5\n",
            "asv_sequences.tsv": "study_id\tasv_id\tsequence\tlength\ttotal_reads\nsynthetic\tA\tACGT\t4\t5\n",
            "asv_sequences.fasta": ">A\nACGT\n",
            "read_tracking.tsv": "study_id\tsample_id\trun_accession\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\tasv_count\nsynthetic\ts1\tr1\t10\t8\t7\t7\t6\t5\t1\n",
            "sample_metadata.tsv": "study_id\tsample_id\trun_accession\nsynthetic\ts1\tr1\n",
            "summary.tsv": "samples\tnonchim_asvs\tmerged_asvs\tinput_pairs\tfiltered_pairs\tmerged_pairs\tnonchim_pairs\n1\t1\t1\t10\t8\t6\t5\n",
            "sequence_length_distribution.tsv": "stage\tlength\tasv_count\tread_count\nmerged\t4\t1\t6\nnonchim\t4\t1\t5\n",
            "error_learning.tsv": "read_direction\tconverged_before_limit\nF\tTRUE\nR\tTRUE\n",
            "input_checksums.tsv": "path\tmd5\ndata/raw/absent.fastq.gz\t" + "f" * 32 + "\n",
            "config_snapshot.yaml": "{}\n", "software_versions.tsv": "software\tversion\nfixture\t1\n",
        }
        for name, text in content.items(): (self.directory / name).write_text(text)
        (self.directory / "SUCCESS").write_text(self.revision + "\n")
        (self.directory / "output_checksums.tsv").write_text("path\tmd5\n" + "".join(f"{name}\t{hashlib.md5((self.directory / name).read_bytes()).hexdigest()}\n" for name in content))
        spec = importlib.util.spec_from_file_location("native_fixture_validation", native_path)
        native_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(native_module)
        native = native_module.validate_run(self.directory)
        native["inputs_checked"] = 99  # Saved historical verification; no raw input access in this review.
        (self.directory / "validation.json").write_text(json.dumps(native))
        cohort = {"status": "valid", "study_id": "synthetic", "samples": 1, "promotion_status": "needs_review"}
        (self.audit / "cohort_validation.json").write_text(json.dumps(cohort))
        (self.audit / "SUCCESS").write_text(self.revision + "\n")
        (self.audit / "validation.json").write_text(json.dumps({**cohort, "source_git_commit": self.revision, "native_validation": native}))
        inputs = [*SCAN.file_records([self.source, native_path], self.root),
                  {"path": "data/raw/absent.fastq.gz", "bytes": 999, "sha256": "f" * 64}]
        (self.audit / "input_checksums.json").write_text(json.dumps(inputs))
        self.outputs = [*self.directory.iterdir(), self.audit / "cohort_validation.json", self.audit / "input_checksums.json"]
        self.record = {"status": "success", "git_commit": self.revision, "parameters": self.policy,
                       "outputs": SCAN.file_records(self.outputs, self.root), "inputs": inputs}
        self.write_record()

    def write_record(self):
        (self.audit / "provenance.json").write_text(json.dumps(self.record))

    def check(self):
        return SCAN.checked_source(self.policy, self.scope, [self.source], self.root)

    def test_valid_native_and_outer_audits_do_not_reopen_absent_fastq(self):
        checked, paths = self.check()
        self.assertEqual(checked["status"], "valid")
        self.assertFalse(checked["fastq_rehashed"])
        self.assertEqual(checked["native_validation"]["inputs_checked"], 0)
        self.assertFalse(any(path.name.endswith(".fastq.gz") for path in paths))

    def test_changed_asv_output_or_source_profile_fails(self):
        for path in (self.directory / "asv_counts.tsv", self.source):
            previous = path.read_bytes()
            path.write_bytes(previous + b"\n")
            with self.subTest(path=path), self.assertRaises(ValueError): self.check()
            path.write_bytes(previous)

    def test_native_validation_must_be_bound_in_outer_outputs(self):
        self.record["outputs"] = [row for row in self.record["outputs"] if not row["path"].endswith("/validation.json")]
        self.write_record()
        with self.assertRaisesRegex(ValueError, "essential DADA2 outputs"): self.check()

    def test_native_validator_must_match_recorded_source_before_import(self):
        path = self.root / "scripts/dada2/validate_outputs.py"
        path.write_text("raise RuntimeError('This changed validator must never be imported')\n")
        with self.assertRaisesRegex(ValueError, "native validator changed"):
            self.check()

    def test_hash_consistent_input_list_must_still_equal_outer_provenance(self):
        changed = deepcopy(self.record["inputs"])
        changed[0]["sha256"] = "0" * 64
        (self.audit / "input_checksums.json").write_text(json.dumps(changed))
        self.record["outputs"] = SCAN.file_records(self.outputs, self.root)
        self.write_record()
        with self.assertRaisesRegex(ValueError, "input checksum list differs"):
            self.check()

    def test_outer_validation_must_match_cohort_and_native_records(self):
        path = self.audit / "validation.json"
        value = json.loads(path.read_text())
        value["samples"] = 2
        path.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, "incomplete or inconsistent"): self.check()


if __name__ == "__main__":
    unittest.main()
