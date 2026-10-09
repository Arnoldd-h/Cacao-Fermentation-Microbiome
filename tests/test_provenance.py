from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.provenance import file_records, write_provenance


class ProvenanceTests(unittest.TestCase):
    def test_records_content_and_dirty_git_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            subprocess.run(["git", "init", "--quiet", str(root)], check=True)
            source = root / "input.txt"
            source.write_bytes(b"observed input")
            subprocess.run(["git", "add", "input.txt"], cwd=root, check=True)
            subprocess.run(["git", "-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "--quiet", "-m", "test fixture"], cwd=root, check=True)
            config = root / "config.json"
            config.write_text("{}", encoding="utf-8")
            output = root / "result.tsv"
            output.write_text("result\n", encoding="utf-8")
            sidecar = root / "result.provenance.json"
            write_provenance(sidecar, root=root, inputs=[source], outputs=[output], config_paths=[config], software_versions={"test": "1"})
            record = json.loads(sidecar.read_text(encoding="utf-8"))
            self.assertTrue(record["git_dirty"])
            self.assertEqual(len(record["git_commit"]), 40)
            self.assertEqual(record["inputs"][0]["path"], "input.txt")
            self.assertEqual(record["inputs"][0]["sha256"], hashlib.sha256(b"observed input").hexdigest())
            self.assertEqual(record["config"][0]["bytes"], 2)
            self.assertEqual(record["outputs"][0]["path"], "result.tsv")
            self.assertIn("+00:00", record["created_at_utc"])

    def test_missing_input_fails_instead_of_inventing_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(FileNotFoundError):
                file_records(["absent"], Path(temporary))


if __name__ == "__main__":
    unittest.main()
