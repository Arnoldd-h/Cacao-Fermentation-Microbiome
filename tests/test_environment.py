from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "environment"))

from verify_environment import parse_r_package_versions  # noqa: E402
from run_in_environment import environment_name, select_prefix  # noqa: E402


class EnvironmentVerificationTests(unittest.TestCase):
    def test_r_package_versions_are_parsed(self) -> None:
        observed = parse_r_package_versions("dada2\t1.34.0\nphyloseq\t1.50.0\n")
        self.assertEqual(observed, [("dada2", "1.34.0"), ("phyloseq", "1.50.0")])

    def test_invalid_r_output_fails(self) -> None:
        with self.assertRaisesRegex(ValueError, "Invalid R package version line"):
            parse_r_package_versions("dada2 1.34.0\n")


class EnvironmentResolutionTests(unittest.TestCase):
    def test_registered_prefix_outside_default_root_is_found(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "alternate" / "cacao-microbiome"
            (prefix / "conda-meta").mkdir(parents=True)
            (prefix / "conda-meta" / "history").touch()
            self.assertEqual(select_prefix("cacao-microbiome", [str(prefix)]), prefix.resolve())

    def test_ambiguous_prefixes_require_explicit_choice(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            prefixes = [Path(directory) / root / "cacao-microbiome" for root in ("a", "b")]
            for prefix in prefixes:
                (prefix / "conda-meta").mkdir(parents=True)
                (prefix / "conda-meta" / "history").touch()
            with self.assertRaisesRegex(ValueError, "Multiple environments"):
                select_prefix("cacao-microbiome", [str(path) for path in prefixes])
            self.assertEqual(select_prefix("cacao-microbiome", [], str(prefixes[1])), prefixes[1].resolve())

    def test_invalid_override_never_falls_back(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "Not an existing"):
                select_prefix("cacao-microbiome", [], directory)

    def test_missing_environment_fails(self) -> None:
        with self.assertRaisesRegex(ValueError, "No registered environment"):
            select_prefix("cacao-microbiome", [])

    def test_name_comes_from_environment_specification(self) -> None:
        self.assertEqual(environment_name(ROOT / "environment" / "environment.yml"), "cacao-microbiome")


if __name__ == "__main__":
    unittest.main()
