"""Synthetic position records distinguish terminal and internal primer matches."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("primer_positions", ROOT / "scripts/study/review_primer_positions.py")
REVIEW = importlib.util.module_from_spec(spec)
spec.loader.exec_module(REVIEW)


class PrimerPositionTests(unittest.TestCase):
    def parse(self, text):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.info.tsv"
            path.write_text(text)
            return REVIEW.read_positions(path)

    def test_terminal_internal_and_unmatched_records(self):
        total, counts = self.parse("a\t0\t2\t4\tAA\tCC\t\tprobe\n"
                                   "b\t1\t2\t4\tAA\tCT\tGG\tprobe\n"
                                   "c\t-1\tACGT\tIIII\n")
        self.assertEqual(total, 3)
        self.assertEqual(dict(counts), {(2, 4, 4): 1, (2, 4, 6): 1})

    def test_invalid_coordinates_and_quality_lengths_fail(self):
        for line in ("a\t0\t3\t5\tAA\tCC\t\tprobe\n", "a\t-1\tAA\tI\n", "a\t0\n"):
            with self.subTest(line=line), self.assertRaises(ValueError):
                self.parse(line)
