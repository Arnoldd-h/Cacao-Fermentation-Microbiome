#!/usr/bin/env python3
"""Verify command-line tools and R packages in the active project environment.

Inputs: executables available on PATH in the active environment.
Output: environment/software_versions.tsv by default.
Fails if a required executable is missing, a version command fails, or an R
package cannot be loaded.
"""

from __future__ import annotations

import argparse
import csv
import os
import platform
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R_PACKAGES = (
    "RcppParallel",
    "dada2",
    "phyloseq",
    "CVXR",
    "Matrix",
    "lme4",
    "reformulas",
    "lmerTest",
    "ANCOMBC",
    "vegan",
    "metafor",
    "yaml",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "environment" / "software_versions.tsv",
    )
    return parser.parse_args()


def run_version(command: list[str]) -> str:
    """Run a version command and return its first non-empty output line."""

    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    combined = "\n".join(part for part in (completed.stdout, completed.stderr) if part)
    for line in combined.splitlines():
        if line.strip():
            return line.strip()
    raise RuntimeError(f"Version command produced no output: {' '.join(command)}")


def parse_r_package_versions(output: str) -> list[tuple[str, str]]:
    """Parse tab-delimited package versions emitted by the R smoke test."""

    versions: list[tuple[str, str]] = []
    for line in output.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 2 or not all(part.strip() for part in parts):
            raise ValueError(f"Invalid R package version line: {line}")
        versions.append((parts[0].strip(), parts[1].strip()))
    return versions


def r_package_versions() -> list[tuple[str, str]]:
    """Load every required R package and return its installed version."""

    package_vector = ",".join(f'"{package}"' for package in R_PACKAGES)
    expression = (
        f"packages <- c({package_vector}); "
        "for (pkg in packages) { "
        "suppressPackageStartupMessages(library(pkg, character.only=TRUE)); "
        "cat(pkg, as.character(packageVersion(pkg)), sep='\\t'); cat('\\n') }"
    )
    try:
        completed = subprocess.run(
            ["Rscript", "-e", expression], check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "R package smoke test failed").strip()
        raise RuntimeError(detail) from exc
    versions = parse_r_package_versions(completed.stdout)
    loaded = {name for name, _ in versions}
    missing = sorted(set(R_PACKAGES) - loaded)
    if missing:
        raise RuntimeError(f"R packages did not report versions: {', '.join(missing)}")
    return versions


def build_rows() -> list[dict[str, str]]:
    """Run all smoke tests and return stable version rows."""

    command_versions = [
        ("Python", platform.python_version(), "runtime loaded"),
        ("FastQC", run_version(["fastqc", "--version"]), "CLI smoke test passed"),
        ("MultiQC", run_version(["multiqc", "--version"]), "CLI smoke test passed"),
        ("Cutadapt", run_version(["cutadapt", "--version"]), "CLI smoke test passed"),
        ("Snakemake", run_version(["snakemake", "--version"]), "CLI smoke test passed"),
        ("R", run_version(["R", "--version"]), "runtime loaded"),
    ]
    rows = [
        {
            "component": component,
            "version": version,
            "validation": validation,
            "platform": "linux-64 (WSL2)",
        }
        for component, version, validation in command_versions
    ]
    rows.extend(
        {
            "component": package,
            "version": version,
            "validation": "R package loaded",
            "platform": "linux-64 (WSL2)",
        }
        for package, version in r_package_versions()
    )
    return rows


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    """Write the version report atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("component", "version", "validation", "platform"),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    args = parse_args()
    rows = build_rows()
    write_rows(args.output, rows)
    print(f"Verified {len(rows)} environment components")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
