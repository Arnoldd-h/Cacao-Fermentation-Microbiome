"""Incomplete public downloads cannot be silently promoted or excluded."""

import gzip
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("audit_downloads", ROOT / "scripts/study/audit_downloads.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class StudyDownloadTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.data = gzip.compress(b"@SRR1.1 1/1\nACGT\n+\nIIII\n", mtime=0)
        self.record = dict(target_relative="fixture/SRR1/SRR1_1.fastq.gz", source_url="https://example.invalid/read.fastq.gz",
                           expected_bytes=str(len(self.data)), expected_md5=hashlib.md5(self.data, usedforsecurity=False).hexdigest())
        self.path = self.root / "data/raw" / self.record["target_relative"]
        self.path.parent.mkdir(parents=True)

    def test_cached_valid_data_need_no_network(self):
        self.path.write_bytes(self.data)
        def unexpected(*args, **kwargs):
            self.fail("Valid raw data must be reused")
        report, paths = module.audit_streams([self.record], self.root, unexpected)
        self.assertEqual(report["status"], "complete")
        self.assertEqual(report["streams"][0]["read_count"], 1)
        self.assertEqual(paths, [self.path])

    def test_matching_md5_does_not_hide_invalid_fastq_structure(self):
        data = gzip.compress(b"@SRR1.1\nACGT\n+\nI\n", mtime=0)
        self.path.write_bytes(data)
        record = {**self.record, "expected_bytes": str(len(data)), "expected_md5": hashlib.md5(data, usedforsecurity=False).hexdigest()}
        report, _ = module.audit_streams([record], self.root)
        self.assertEqual(report["streams"][0]["status"], "invalid_local_file")

    def test_html_response_blocks_without_exclusion_or_file_creation(self):
        class Response:
            status = 200
            headers = {"Content-Type": "text/html;charset=UTF-8"}
            def geturl(self): return "https://example.invalid/read.fastq.gz/"
            def __enter__(self): return self
            def __exit__(self, *args): return False
        report, paths = module.audit_streams([self.record], self.root, lambda *a, **kw: Response())
        self.assertEqual(report["status"], "blocked")
        self.assertEqual(report["streams"][0]["status"], "html_instead_of_fastq")
        self.assertEqual(report["processing_exclusions_added"], 0)
        self.assertFalse(self.path.exists())
        self.assertEqual(paths, [])

    def test_altered_raw_is_preserved_and_reported(self):
        self.path.write_bytes(b"altered")
        report, _ = module.audit_streams([self.record], self.root)
        self.assertEqual(report["streams"][0]["status"], "invalid_local_file")
        self.assertEqual(self.path.read_bytes(), b"altered")


if __name__ == "__main__":
    unittest.main()
