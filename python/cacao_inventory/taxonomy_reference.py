"""Pinned DADA2 reference validation and streaming FASTA inspection."""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from urllib.parse import urlparse


def repository_path(root: Path, label: str) -> Path:
    path = (root / label).resolve()
    if Path(label).is_absolute() or root.resolve() not in path.parents:
        raise ValueError("Reference path must remain inside the repository")
    return path


def validate_taxonomy_config(config: dict, root: Path) -> None:
    reference = config["reference"]
    path = repository_path(root, reference["path"])
    if not path.is_relative_to((root / "references/databases").resolve()):
        raise ValueError("Reference database must remain under references/databases")
    repository_path(root, reference["provenance_path"])
    if urlparse(reference["url"]).scheme != "https":
        raise ValueError("Reference URL must use HTTPS")
    if not re.fullmatch(r"[a-f0-9]{32}", reference["expected_md5"]):
        raise ValueError("Invalid reference MD5")
    if type(reference["expected_bytes"]) is not int or reference["expected_bytes"] <= 0:
        raise ValueError("Invalid reference byte count")
    settings = config["classification"]
    if settings.get("raw_min_boot") != 0 or settings.get("rng_kind") != ["Mersenne-Twister", "Inversion", "Rejection"]:
        raise ValueError("Unfiltered calls and the preregistered RNG must be explicit")
    if settings["tax_levels"] != ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]:
        raise ValueError("Expected the six original SILVA ranks")
    thresholds = [settings["primary_min_boot"], *settings["sensitivity_min_boot"]]
    if len(set(thresholds)) != len(thresholds) or any(type(t) is not int or not 1 <= t <= 100 for t in thresholds):
        raise ValueError("Invalid or duplicate bootstrap thresholds")
    if settings["method"] != "dada2::assignTaxonomy" or settings["species_assignment"] is not False:
        raise ValueError("Only genus-level assignTaxonomy is implemented")
    if settings["output_bootstraps"] is not True or type(settings["try_reverse_complement"]) is not bool:
        raise ValueError("Bootstrap output is required; reverse-complement policy must be explicit")
    for field in ("random_seed", "threads"):
        if type(settings[field]) is not int or settings[field] < (1 if field == "threads" else 0):
            raise ValueError(f"Invalid {field}")
    if config["screening"]["exclude_flagged_asvs"] is not False:
        raise ValueError("Exclusion of flagged ASVs is not implemented")
    re.compile(config["screening"]["unresolved_genus_pattern"])
    download = config["download"]
    for field in ("timeout_seconds", "retries", "chunk_size_bytes", "minimum_free_space_factor"):
        value = download[field]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0 or (
                field != "minimum_free_space_factor" and type(value) is not int):
            raise ValueError(f"Invalid download {field}")


def inspect_training_fasta(path: Path, levels: int) -> dict[str, int]:
    """Check gzip integrity, nonempty records and lineage depth without loading it."""
    records = bases = sequence_bases = 0
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if records and not sequence_bases:
                    raise ValueError("Empty reference FASTA sequence")
                lineage = line[1:].split(";")
                if lineage[-1] == "":
                    lineage.pop()
                if not 1 <= len(lineage) <= levels or not lineage[0]:
                    raise ValueError("Reference lineage has incompatible ranks")
                records += 1
                sequence_bases = 0
            else:
                if not records or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", line):
                    raise ValueError("Invalid reference FASTA sequence")
                bases += len(line)
                sequence_bases += len(line)
    if not records or not sequence_bases:
        raise ValueError("Empty reference FASTA")
    return {"reference_sequences": records, "reference_bases": bases}
