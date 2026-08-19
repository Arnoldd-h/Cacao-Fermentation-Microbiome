from __future__ import annotations

import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.download import download_file, expand_pilot_manifest  # noqa: E402
from cacao_inventory.schema import PILOT_MANIFEST_COLUMNS  # noqa: E402


def manifest_row() -> dict[str, str]:
    row = {column: "" for column in PILOT_MANIFEST_COLUMNS}
    row.update(
        {
            "study_id": "pilot_study",
            "bioproject": "PRJNA123456",
            "sample_id": "sample_1",
            "run_accession": "SRR1234567",
            "fastq_ftp": "https://example.org/read_1.fastq.gz;https://example.org/read_2.fastq.gz",
            "fastq_md5": "a" * 32 + ";" + "b" * 32,
            "fastq_bytes": "10;20",
        }
    )
    return row


class DownloadTests(unittest.TestCase):
    def test_manifest_expands_to_unique_fastq_records(self) -> None:
        records = expand_pilot_manifest([manifest_row()])

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["read_file"], "read_1.fastq.gz")
        self.assertEqual(records[1]["expected_bytes"], "20")
        self.assertNotEqual(records[0]["target_relative"], records[1]["target_relative"])

    def test_unaligned_manifest_fields_fail(self) -> None:
        row = manifest_row()
        row["fastq_md5"] = "a" * 32

        with self.assertRaisesRegex(ValueError, "Unaligned FASTQ metadata"):
            expand_pilot_manifest([row])

    def test_valid_existing_file_is_reused_without_network(self) -> None:
        content = b"validated fastq fixture"
        checksum = hashlib.md5(content, usedforsecurity=False).hexdigest()
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "fixture.fastq.gz"
            destination.write_bytes(content)

            status = download_file(
                "https://example.invalid/fixture.fastq.gz",
                destination,
                len(content),
                checksum,
                1024,
            )

        self.assertEqual(status, "reused_valid")

    def test_invalid_existing_file_is_not_overwritten(self) -> None:
        content = b"user data remains intact"
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "fixture.fastq.gz"
            destination.write_bytes(content)

            with self.assertRaisesRegex(ValueError, "MD5 mismatch"):
                download_file(
                    "https://example.invalid/fixture.fastq.gz",
                    destination,
                    len(content),
                    "0" * 32,
                    1024,
                )
            self.assertEqual(destination.read_bytes(), content)


if __name__ == "__main__":
    unittest.main()
