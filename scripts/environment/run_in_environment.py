#!/usr/bin/env python3
"""Run a command in the registered project Micromamba environment.

Inputs: environment/environment.yml (name), Micromamba's registered prefixes,
and an optional --prefix or CACAO_ENV_PREFIX override. Outputs are those of
the command; no environment is installed or modified. Fails if the prefix is
missing or ambiguous, Micromamba is unavailable, or the command fails.
Run inside Linux/WSL, for example: python3 scripts/environment/run_in_environment.py
snakemake --snakefile workflow/Snakefile --cores 2.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOGGER = logging.getLogger(__name__)


def environment_name(specification: Path) -> str:
    """Read the single top-level name without requiring a YAML dependency."""
    matches = re.findall(r"^name:\s*([A-Za-z0-9_.-]+)\s*$", specification.read_text(), re.M)
    if len(matches) != 1:
        raise ValueError(f"Expected one simple environment name in {specification}")
    return matches[0]


def select_prefix(name: str, registered: list[str], override: str | None = None) -> Path:
    """Resolve a unique existing environment, rejecting silent fallback."""
    if override:
        prefix = Path(override).expanduser().resolve()
        if not (prefix / "conda-meta" / "history").is_file():
            raise ValueError(f"Not an existing Conda environment: {prefix}")
        return prefix
    matches = {
        Path(value).expanduser().resolve()
        for value in registered
        if Path(value).name == name
        and (Path(value).expanduser() / "conda-meta" / "history").is_file()
    }
    if not matches:
        raise ValueError(
            f"No registered environment named {name!r}; create it from the lock file "
            "or set CACAO_ENV_PREFIX to its existing prefix"
        )
    if len(matches) > 1:
        raise ValueError(
            f"Multiple environments named {name!r}; choose one with --prefix "
            "or CACAO_ENV_PREFIX"
        )
    return next(iter(matches))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", help="Explicit existing Linux environment prefix")
    parser.add_argument("command", nargs=argparse.REMAINDER, help="Command and its arguments")
    args = parser.parse_args()
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("a command is required")
    return args


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    executable = shutil.which("micromamba")
    if executable is None:
        raise RuntimeError("micromamba is unavailable on PATH; run inside the configured WSL shell")
    override = args.prefix or os.environ.get("CACAO_ENV_PREFIX")
    registered: list[str] = []
    if not override:
        completed = subprocess.run(
            [executable, "env", "list", "--json"], check=True, capture_output=True, text=True
        )
        payload = json.loads(completed.stdout)
        registered = payload.get("envs", [])
        if not isinstance(registered, list) or not all(isinstance(item, str) for item in registered):
            raise ValueError("Unexpected Micromamba environment-list schema")
    prefix = select_prefix(environment_name(ROOT / "environment" / "environment.yml"), registered, override)
    LOGGER.info("Using environment: %s", prefix)
    completed = subprocess.run(
        [executable, "run", "--prefix", str(prefix), *args.command], cwd=ROOT
    )
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
