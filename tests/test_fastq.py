from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.fastq import validate_fastq_file  # noqa: E402


class FastqValidationTests(unittest.TestCase):
    def test_valid_gzip_fastq_reports_counts(self) -> None:
        content = b"@read1\nACGT\n+\nIIII\n@read2\nAC\n+\nII\n"
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fixture.fastq.gz"
            with gzip.open(path, "wb") as handle:
                handle.write(content)
            metrics = validate_fastq_file(path)

        self.assertEqual(metrics["read_count"], 2)
        self.assertEqual(metrics["base_count"], 6)
        self.assertEqual(metrics["minimum_read_length"], 2)
        self.assertEqual(metrics["maximum_read_length"], 4)

    def test_sequence_quality_mismatch_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fixture.fastq.gz"
            with gzip.open(path, "wb") as handle:
                handle.write(b"@read1\nACGT\n+\nIII\n")

            with self.assertRaisesRegex(ValueError, "Sequence/quality length mismatch"):
                validate_fastq_file(path)

    def test_truncated_record_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fixture.fastq.gz"
            with gzip.open(path, "wb") as handle:
                handle.write(b"@read1\nACGT\n+\n")

            with self.assertRaisesRegex(ValueError, "Truncated FASTQ record"):
                validate_fastq_file(path)


if __name__ == "__main__":
    unittest.main()
