#!/usr/bin/env python3
"""Measure verified primer constructs in every pilot FASTQ with Cutadapt.

Inputs: project configuration and the six-run pilot manifest. The raw FASTQ
paths are resolved below data/raw from study and run identifiers. Output:
results/qc/pilot/primer_detection.tsv by default. Fails on invalid IUPAC,
missing paired FASTQ, Cutadapt errors, or inconsistent read counts.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from cacao_inventory.config import load_json_yaml  # noqa: E402
from cacao_inventory.io import read_tsv, write_tsv_atomic  # noqa: E402
from cacao_inventory.qc import (  # noqa: E402
    cutadapt_detection_counts,
    reverse_complement_iupac,
    validate_iupac_sequence,
)
from cacao_inventory.schema import (  # noqa: E402
    PILOT_MANIFEST_COLUMNS,
    PRIMER_DETECTION_COLUMNS,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "config.yaml")
    parser.add_argument(
        "--manifest", type=Path, default=ROOT / "metadata" / "pilot_manifest.tsv"
    )
    parser.add_argument("--input-stage", choices=("raw", "trimmed"), default="raw")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cutadapt", default="cutadapt")
    return parser.parse_args()


def _processing_config(config: dict[str, Any], bioproject: str) -> dict[str, Any]:
    try:
        processing = config["amplicon_processing"][bioproject]
        primers = processing["primers"]
        detection = processing["primer_detection"]
    except (KeyError, TypeError) as exc:
        raise ValueError(f"Missing amplicon primer configuration for {bioproject}") from exc
    if set(primers) != {"forward", "reverse"}:
        raise ValueError(f"{bioproject} must configure forward and reverse primers")
    try:
        error_rate = float(detection["error_rate"])
        overlap_fraction = float(detection["minimum_overlap_fraction"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid primer detection parameters for {bioproject}") from exc
    if not 0 <= error_rate <= 1 or overlap_fraction != 1.0:
        raise ValueError(
            "Primer detection requires error_rate within [0,1] and full-length overlap"
        )
    if detection.get("allow_indels") is not False:
        raise ValueError("Primer detection must explicitly disable indels")
    return processing


def _run_cutadapt(
    executable: str,
    read_file: Path,
    adapter_name: str,
    search_sequence: str,
    anchored: bool,
    error_rate: float,
) -> tuple[int, int]:
    adapter = f"{adapter_name}={'^' if anchored else ''}{search_sequence}"
    with tempfile.TemporaryDirectory(prefix="cacao-primer-") as temporary:
        report_path = Path(temporary) / "cutadapt.json"
        command = [
            executable,
            "--action=none",
            "--no-indels",
            "--error-rate",
            str(error_rate),
            "--overlap",
            str(len(search_sequence)),
            "--front",
            adapter,
            "--json",
            str(report_path),
            "--output",
            os.devnull,
            str(read_file),
        ]
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            message = completed.stderr.strip() or completed.stdout.strip()
            raise RuntimeError(f"Cutadapt failed for {read_file.name}: {message}")
        with report_path.open(encoding="utf-8") as handle:
            return cutadapt_detection_counts(json.load(handle))


def main() -> int:
    args = parse_args()
    output = args.output or (
        ROOT
        / "results"
        / "qc"
        / "pilot"
        / (
            "primer_detection.tsv"
            if args.input_stage == "raw"
            else "primer_detection_trimmed.tsv"
        )
    )
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    config = load_json_yaml(args.config)
    manifest = read_tsv(args.manifest, PILOT_MANIFEST_COLUMNS)
    bioprojects = {row["bioproject"] for row in manifest}
    if len(bioprojects) != 1:
        raise ValueError("Pilot primer detection requires exactly one BioProject")
    bioproject = next(iter(bioprojects))
    processing = _processing_config(config, bioproject)
    detection = processing["primer_detection"]
    error_rate = float(detection["error_rate"])
    rows: list[dict[str, str]] = []

    for sample in manifest:
        for direction, suffix in (("R1", "1"), ("R2", "2")):
            if args.input_stage == "raw":
                read_file = (
                    ROOT
                    / "data"
                    / "raw"
                    / sample["study_id"]
                    / sample["run_accession"]
                    / f"{sample['run_accession']}_{suffix}.fastq.gz"
                )
            else:
                read_file = (
                    ROOT
                    / "data"
                    / "interim"
                    / sample["study_id"]
                    / "pilot"
                    / sample["run_accession"]
                    / f"{sample['run_accession']}_{suffix}.fastq.gz"
                )
            if not read_file.is_file():
                raise FileNotFoundError(f"Missing pilot FASTQ: {read_file}")
            for primer_role, primer in processing["primers"].items():
                primer_label = "515F" if primer_role == "forward" else "806R"
                for sequence_type, config_key, anchored in (
                    ("gene_specific", "gene_specific_sequence", False),
                    ("full_construct", "full_construct_sequence", True),
                ):
                    primer_sequence = validate_iupac_sequence(
                        str(primer[config_key]),
                        f"{bioproject}.{primer_role}.{config_key}",
                    )
                    for orientation, search_sequence in (
                        ("as_synthesized", primer_sequence),
                        ("reverse_complement", reverse_complement_iupac(primer_sequence)),
                    ):
                        expected = (
                            orientation == "as_synthesized"
                            and (
                                (primer_role == "forward" and direction == "R1")
                                or (primer_role == "reverse" and direction == "R2")
                            )
                        )
                        examined, matched = _run_cutadapt(
                            args.cutadapt,
                            read_file,
                            f"{primer_label}_{sequence_type}_{orientation}",
                            search_sequence,
                            anchored,
                            error_rate,
                        )
                        rows.append(
                            {
                                "study_id": sample["study_id"],
                                "bioproject": bioproject,
                                "sample_id": sample["sample_id"],
                                "run_accession": sample["run_accession"],
                                "read_direction": direction,
                                "read_file": read_file.name,
                                "primer": primer_label,
                                "sequence_type": sequence_type,
                                "primer_sequence": primer_sequence,
                                "search_sequence": search_sequence,
                                "reads_examined": str(examined),
                                "primer_matches": str(matched),
                                "match_percent": f"{matched / examined * 100:.6f}",
                                "orientation": orientation,
                                "search_scope": (
                                    "five_prime_anchored"
                                    if anchored
                                    else "any_position_full_length"
                                ),
                                "error_rate": str(error_rate),
                                "minimum_overlap": str(len(search_sequence)),
                                "expected_case": str(expected).lower(),
                            }
                        )

    write_tsv_atomic(output, rows, PRIMER_DETECTION_COLUMNS)
    logging.info(
        "Recorded %d primer detection tests across %d FASTQ files",
        len(rows),
        len(manifest) * 2,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
