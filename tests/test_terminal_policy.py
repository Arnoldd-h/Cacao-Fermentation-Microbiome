"""Registered terminal method and real Cutadapt behavior on synthetic reads."""
import copy
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/terminal_qc"))
import run as terminal


class TerminalPolicyTests(unittest.TestCase):
    def test_registered_policy_and_separate_paths(self):
        policy, scope, manifest, processing = terminal.load_policy(ROOT / "config/prjna627078_terminal_qc.yaml")
        self.assertEqual(len(manifest), scope["expected_runs"])
        self.assertEqual(list(map(len, terminal.sequences(processing))), [21, 17])
        self.assertEqual(processing["dada2"]["filter"]["min_len"], 150)
        inputs, outputs = terminal.read_paths(policy, scope, manifest[0])
        self.assertFalse(set(inputs) & set(outputs))

    def test_policy_rejects_filtering_partial_internal_and_overwrite(self):
        policy = json.loads((ROOT / "config/prjna627078_terminal_qc.yaml").read_text())
        for field, value in (("discard_reads", True), ("partial_matches", True), ("internal_matches", True),
                             ("allow_indels", True), ("passes", 2), ("error_rate", .2),
                             ("interim_dir", "data/interim/mota_gutierrez_2021_mexico/full"),
                             ("interim_dir", "data/interim/another_study/terminal_3prime")):
            altered = copy.deepcopy(policy)
            altered[field] = value
            with self.subTest(field=field), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "policy.json"
                path.write_text(json.dumps(altered))
                with self.assertRaises(ValueError):
                    terminal.load_policy(path)

    def test_run_evidence_rejects_stale_method_missing_and_changed_reports(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            policy = {"cutadapt_version": "5.0"}
            method = [{"path": "synthetic.py", "sha256": "synthetic"}]
            for name in terminal.REPORTS:
                (path / name).write_text("{}")
            audit = {"method_checksums": method, "cutadapt_version": "5.0",
                     "report_checksums": terminal.file_records([path / name for name in terminal.REPORTS], ROOT)}
            with self.assertRaisesRegex(ValueError, "different method"):
                terminal.validate_run_record(audit, path, [], policy)
            (path / terminal.REPORTS[0]).write_text('{"changed": true}')
            with self.assertRaisesRegex(ValueError, "changed or disappeared"):
                terminal.validate_run_record(audit, path, method, policy)
            (path / terminal.REPORTS[0]).unlink()
            with self.assertRaises(FileNotFoundError):
                terminal.validate_run_record(audit, path, method, policy)

    @unittest.skipUnless(shutil.which("cutadapt"), "Cutadapt integration requires the registered environment")
    def test_real_cutadapt_full_suffix_degeneracy_errors_and_pair_conservation(self):
        # R ambiguity permits A/G. N in a read is not treated as an adapter wildcard.
        primer = "ACGTRACGTACGTACGT"
        exact = primer.replace("R", "A")
        one_error = "T" + exact[1:]
        two_errors = "TT" + exact[2:]
        sequences = ["GGGG" + exact, "GGGG" + exact + "CC", "GGGG" + exact[:-1],
                     "GGGG" + primer.replace("R", "G"), "GGGG" + one_error,
                     "GGGG" + two_errors, exact, "GGGG", "GGGG" + exact[:5] + "N" + exact[6:]]
        expected = ["GGGG", sequences[1], sequences[2], "GGGG", "GGGG", sequences[5], "", "GGGG", "GGGG"]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = [root / f"synthetic_{mate}.fastq.gz" for mate in (1, 2)]
            outputs = [root / f"clean_{mate}.fastq.gz" for mate in (1, 2)]
            for mate, path in enumerate(inputs, 1):
                with gzip.open(path, "wt") as handle:
                    for i, sequence in enumerate(sequences):
                        handle.write(f"@synthetic{i}/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            command = terminal.cutadapt_command({"error_rate": .1}, [primer, primer], inputs, outputs, root / "cutadapt.json")
            subprocess.run(command, check=True, capture_output=True, text=True)
            for path in outputs:
                with gzip.open(path, "rt") as handle:
                    actual = handle.read().splitlines()[1::4]
                self.assertEqual(actual, expected)
            audit = terminal.audit_pair(inputs, outputs, [len(primer)] * 2, 4)
            self.assertEqual(audit["pairs"], len(sequences))
            self.assertEqual(audit["directions"]["R1"]["changed_reads"], 5)
            self.assertEqual(audit["pairs_below_min_len"], 1)


if __name__ == "__main__":
    unittest.main()
