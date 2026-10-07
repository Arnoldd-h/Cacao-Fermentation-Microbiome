"""Atomic metadata writes preserve workflow inputs when content is unchanged."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.io import read_tsv, write_tsv_atomic


class AtomicTsvTests(unittest.TestCase):
    def test_identical_content_preserves_mtime_and_removes_temporary_file(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "metadata.tsv"
            rows = [{"sample_id": "sample_1", "country": "México"}]
            columns = ["sample_id", "country"]
            write_tsv_atomic(path, rows, columns)
            os.utime(path, ns=(1_600_000_000_000_000_000, 1_600_000_000_000_000_000))
            before = path.stat().st_mtime_ns
            write_tsv_atomic(path, rows, columns)
            self.assertEqual(path.stat().st_mtime_ns, before)
            self.assertEqual(read_tsv(path, columns), rows)
            self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_changed_content_replaces_existing_file(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "metadata.tsv"
            write_tsv_atomic(path, [{"value": "old"}], ["value"])
            write_tsv_atomic(path, [{"value": "new"}], ["value"])
            self.assertEqual(read_tsv(path, ["value"]), [{"value": "new"}])
            self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_invalid_row_leaves_existing_destination_intact(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "metadata.tsv"
            write_tsv_atomic(path, [{"value": "valid"}], ["value"])
            with self.assertRaises(ValueError):
                write_tsv_atomic(path, [{"extra": "invalid"}], ["value"])
            self.assertEqual(read_tsv(path, ["value"]), [{"value": "valid"}])
            self.assertEqual(list(Path(directory).iterdir()), [path])
