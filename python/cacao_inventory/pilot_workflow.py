"""Pure helpers for the manifest-driven sequencing workflow."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from .download import expand_pilot_manifest
from .io import read_tsv
from .schema import PILOT_MANIFEST_COLUMNS


def load_pilot_manifest(path: str | Path) -> list[dict[str, str]]:
    """Validate the single-study paired pilot without guessing filenames."""

    rows = read_tsv(path, PILOT_MANIFEST_COLUMNS)
    if not rows:
        raise ValueError("Pilot manifest must not be empty")
    for field in ("study_id", "bioproject"):
        if len({row[field] for row in rows}) != 1:
            raise ValueError(f"Pilot manifest requires exactly one {field}")
    if len({row["run_accession"] for row in rows}) != len(rows):
        raise ValueError("Duplicate run_accession in pilot manifest")
    for row in rows:
        for field in ("study_id", "run_accession"):
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_\-]*", row[field]):
                raise ValueError(f"Unsafe pilot {field}: {row[field]}")
        if row["library_layout"] != "PAIRED":
            raise ValueError("Pilot sequencing workflow requires paired FASTQ")
        records = expand_pilot_manifest([row])
        expected = {f"{row['run_accession']}_{direction}.fastq.gz" for direction in (1, 2)}
        if {record["read_file"] for record in records} != expected or len(records) != 2:
            raise ValueError(f"Pilot requires ENA _1/_2 FASTQ names for {row['run_accession']}")
    return rows


def manifest_row(rows: list[dict[str, str]], run_accession: str, study_id: str | None = None) -> dict[str, str]:
    """Find exactly one run, optionally checking the study wildcard."""

    selected = [row for row in rows if row["run_accession"] == run_accession]
    if len(selected) != 1:
        raise ValueError(f"Run is absent or ambiguous in pilot manifest: {run_accession}")
    row = selected[0]
    if study_id is not None and row["study_id"] != study_id:
        raise ValueError(f"Study wildcard disagrees with manifest for {run_accession}")
    return row


def fastq_path(row: dict[str, str], direction: str, stage: str = "raw") -> str:
    """Resolve a raw or interim FASTQ beneath its manifest study and run."""

    if direction not in ("1", "2"):
        raise ValueError(f"Invalid paired direction: {direction}")
    if stage == "raw":
        directory = f"data/raw/{row['study_id']}/{row['run_accession']}"
    elif stage == "trimmed":
        directory = f"data/interim/{row['study_id']}/pilot/{row['run_accession']}"
    else:
        raise ValueError(f"Unknown FASTQ stage: {stage}")
    return f"{directory}/{row['run_accession']}_{direction}.fastq.gz"


def processing_configuration(config: dict[str, Any], row: dict[str, str]) -> dict[str, Any]:
    """Use the BioProject actually selected by the manifest."""

    try:
        processing = config["amplicon_processing"][row["bioproject"]]
    except KeyError as exc:
        raise ValueError(f"Missing amplicon_processing for {row['bioproject']}") from exc
    if processing.get("study_id") != row["study_id"]:
        raise ValueError("Amplicon configuration study_id disagrees with pilot manifest")
    return processing
