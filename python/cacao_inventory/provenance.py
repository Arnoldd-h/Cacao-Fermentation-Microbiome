"""Small, hash-based execution records for reproducible scientific outputs."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def sha256_file(path: str | Path) -> str:
    """Hash a file incrementally without loading sequencing data into memory."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_records(paths: Iterable[str | Path], root: Path) -> list[dict[str, Any]]:
    """Record stable relative paths, byte sizes and SHA-256 digests."""

    records = []
    for value in sorted({str(path) for path in paths}):
        path = Path(value)
        path = (root / path).resolve() if not path.is_absolute() else path.resolve()
        try:
            label = path.relative_to(root.resolve()).as_posix()
        except ValueError:
            label = path.as_posix()
        records.append({"path": label, "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    return records


def collect_versions(*executables: str) -> dict[str, str]:
    """Read actual versions; a missing/failing tool is a provenance error."""

    versions = {"Python": platform.python_version()}
    for executable in executables:
        result = subprocess.run(
            [executable, "--version"], check=True, capture_output=True, text=True, timeout=60
        )
        lines = (result.stdout + "\n" + result.stderr).splitlines()
        versions[executable] = next((line.strip() for line in lines if line.strip()), "")
        if not versions[executable]:
            raise ValueError(f"No version reported by {executable}")
    return versions


def write_provenance(
    path: str | Path,
    *,
    root: str | Path,
    inputs: Iterable[str | Path],
    outputs: Iterable[str | Path],
    config_paths: Iterable[str | Path],
    software_versions: dict[str, str],
    command: list[str] | None = None,
    parameters: dict[str, Any] | None = None,
) -> None:
    """Write an atomic sidecar including Git dirty state and all content hashes."""

    root = Path(root).resolve()
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, check=True, capture_output=True, text=True
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=normal"],
        cwd=root, check=True, capture_output=True, text=True,
    ).stdout
    record = {
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": commit,
        "git_dirty": bool(status.strip()),
        "config": file_records(config_paths, root),
        "inputs": file_records(inputs, root),
        "outputs": file_records(outputs, root),
        "software_versions": software_versions,
        "command": command or [],
        "parameters": parameters or {},
    }
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(dir=destination.parent, suffix=".tmp")
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, destination)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def qc_provenance(
    output: Path, *, root: Path, inputs: Iterable[str | Path],
    outputs: Iterable[str | Path] | None = None, tools: tuple[str, ...] = (),
    config: Path | None = None, command: list[str] | None = None,
) -> None:
    """Record one QC stage, including its executable source and workflow."""

    code = list((root / "python" / "cacao_inventory").glob("*.py"))
    code += list((root / "scripts" / "qc").glob("*.py"))
    code += [root / "workflow" / "Snakefile", root / "workflow" / "rules" / "pilot_qc.smk"]
    code += [root / "environment" / "conda-linux-64.lock"]
    write_provenance(
        output.with_suffix(".provenance.json"), root=root,
        inputs=[*inputs, *code], outputs=list(outputs) if outputs is not None else [output],
        config_paths=[config or root / "config" / "config.yaml"],
        software_versions=collect_versions(*tools), command=command,
    )
