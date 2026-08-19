"""Atomic TSV input/output utilities."""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path

from .schema import TABLE_SCHEMAS


def write_tsv_atomic(path: str | Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    """Write a TSV atomically so interrupted runs do not leave partial metadata."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp", text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, destination)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def write_inventory(output_dir: str | Path, tables: dict[str, list[dict[str, str]]]) -> None:
    """Write all inventory tables using their stable schemas."""

    target = Path(output_dir)
    for name, columns in TABLE_SCHEMAS.items():
        write_tsv_atomic(target / f"{name}.tsv", tables[name], columns)


def read_tsv(path: str | Path, required_columns: list[str]) -> list[dict[str, str]]:
    """Read a TSV and fail if its header no longer satisfies the schema."""

    source = Path(path)
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = set(reader.fieldnames or [])
        missing = sorted(set(required_columns) - headers)
        if missing:
            raise ValueError(f"{source} is missing required columns: {', '.join(missing)}")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def read_inventory(input_dir: str | Path) -> dict[str, list[dict[str, str]]]:
    """Read all versioned inventory tables."""

    source = Path(input_dir)
    return {
        name: read_tsv(source / f"{name}.tsv", columns)
        for name, columns in TABLE_SCHEMAS.items()
    }
