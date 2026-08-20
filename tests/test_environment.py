from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "environment"))

from verify_environment import parse_r_package_versions  # noqa: E402


class EnvironmentVerificationTests(unittest.TestCase):
    def test_r_package_versions_are_parsed(self) -> None:
        observed = parse_r_package_versions("dada2\t1.34.0\nphyloseq\t1.50.0\n")
        self.assertEqual(observed, [("dada2", "1.34.0"), ("phyloseq", "1.50.0")])

    def test_invalid_r_output_fails(self) -> None:
        with self.assertRaisesRegex(ValueError, "Invalid R package version line"):
            parse_r_package_versions("dada2 1.34.0\n")


if __name__ == "__main__":
    unittest.main()
