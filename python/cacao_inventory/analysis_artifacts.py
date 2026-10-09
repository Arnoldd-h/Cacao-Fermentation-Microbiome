"""Strict matrix input and hash records for small downstream analysis stages."""

from __future__ import annotations

import csv
import json
import math
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from .provenance import file_records
from .taxonomy_reference import repository_path


def read_table(path: Path, required: list[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if not header or len(set(header)) != len(header) or not set(required).issubset(header):
            raise ValueError(f"Invalid or duplicate columns in {path}")
        rows = []
        for values in reader:
            if len(values) != len(header):
                raise ValueError(f"Malformed row in {path}")
            rows.append(dict(zip(header, values)))
        return rows


def unique_rows(rows: list[dict], fields: tuple[str, ...]) -> dict:
    result = {tuple(row[field] for field in fields): row for row in rows}
    if len(result) != len(rows) or any(not all(key) for key in result):
        raise ValueError(f"Empty or duplicate identifiers: {fields}")
    return result


def read_counts(path: Path) -> tuple[list[dict], list[str]]:
    rows = read_table(path, ["study_id", "sample_id"])
    unique_rows(rows, ("study_id", "sample_id"))
    if not rows or len({row["study_id"] for row in rows}) != 1:
        raise ValueError("Counts require one nonempty study")
    features = [field for field in rows[0] if field not in ("study_id", "sample_id")]
    if not features or any(not feature for feature in features):
        raise ValueError("Empty feature columns")
    for row in rows:
        for feature in features:
            if not re.fullmatch(r"[0-9]+", row[feature]):
                raise ValueError("Counts must be nonnegative integers")
    if any(sum(int(row[feature]) for row in rows) == 0 for feature in features):
        raise ValueError("All-zero feature in count matrix")
    return rows, features


def validate_policy(config: dict, root: Path) -> dict[str, Path]:
    paths = {key: repository_path(root, value) for key, value in config["paths"].items()}
    if set(paths) != {"dada2_dir", "taxonomy_dir", "bacterial_dir", "filtering_dir", "diversity_dir"} or len(set(paths.values())) != 5:
        raise ValueError("Stage paths must be distinct and complete")
    if not paths["bacterial_dir"].is_relative_to(root.resolve() / "data/processed"):
        raise ValueError("Derived counts must remain under data/processed")
    for key in ("filtering_dir", "diversity_dir"):
        if not paths[key].is_relative_to(root.resolve() / "results"):
            raise ValueError("Reports must remain under results")
    if any(a in b.parents for a in paths.values() for b in paths.values() if a != b):
        raise ValueError("Stage directories cannot contain one another")
    filt, analysis = config["filtering"], config["analysis"]
    if filt["taxonomy_call"] != "primary" or filt["tax_levels"] != ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]:
        raise ValueError("Only primary six-rank taxonomy is supported")
    if type(filt["minimum_bootstrap"]) is not int or not 1 <= filt["minimum_bootstrap"] <= 100:
        raise ValueError("Invalid taxonomy support threshold")
    if filt["unclassified_kingdom_policy"] != "retain_and_flag" or filt["empty_sample_policy"] != "fail" or filt["low_abundance_filter"] is not False:
        raise ValueError("Unsupported exclusion or empty-sample policy")
    if type(filt["exclude_known_non_bacterial"]) is not bool or not filt["bacterial_kingdom"]:
        raise ValueError("Invalid Kingdom policy")
    labels = filt["organelle_labels"]
    if not labels or any(not isinstance(label, str) or not label for label in labels) or len({label.casefold() for label in labels}) != len(labels):
        raise ValueError("Invalid organelle labels")
    expected = {"feature_level": "within_study_asv", "sample_aggregation": "none", "inference": False,
                "rarefaction": False, "log_base": "natural", "primary_beta": "euclidean_clr",
                "zero_policy": "add_pseudocount_to_all_counts", "secondary_beta": "bray_curtis_relative_abundance",
                "pca_center": True, "pca_scale": False,
                "alpha_metrics": ["observed_asvs", "shannon", "gini_simpson", "inverse_simpson", "shannon_effective"]}
    if any(analysis[key] != value or (isinstance(value, bool) and type(analysis[key]) is not bool) for key, value in expected.items()):
        raise ValueError("Unsupported diversity method or inference requested")
    constants = [analysis["primary_pseudocount"], *analysis["sensitivity_pseudocounts"]]
    if not constants or any(isinstance(c, bool) or not isinstance(c, (int, float)) or not math.isfinite(c) or c <= 0 for c in constants) or len(set(constants)) != len(constants):
        raise ValueError("Pseudocounts must be finite, positive and unique")
    if type(analysis["random_seed"]) is not int or analysis["random_seed"] < 0:
        raise ValueError("Invalid seed")
    tolerance = analysis["validation_tolerance"]
    if isinstance(tolerance, bool) or not isinstance(tolerance, (int, float)) or not 0 < tolerance <= 1e-6:
        raise ValueError("Invalid validation tolerance")
    figures = config["figures"]
    if type(figures["png_dpi"]) is not int or figures["png_dpi"] < 300:
        raise ValueError("Publication PNG requires at least 300 dpi")
    for key in ("width_inches", "height_inches"):
        if isinstance(figures[key], bool) or not isinstance(figures[key], (int, float)) or not math.isfinite(figures[key]) or figures[key] <= 0:
            raise ValueError("Invalid figure dimensions")
    if any(not re.fullmatch(r"#[0-9a-fA-F]{6}", color) for color in figures["stage_colors"].values()):
        raise ValueError("Invalid stage color")
    return paths


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def begin_stage(report: Path, root: Path) -> dict:
    report.mkdir(parents=True, exist_ok=True)
    state = {"git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
             "git_status_at_start": subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).splitlines(),
             "started_at_utc": datetime.now(timezone.utc).isoformat()}
    state["git_dirty"] = bool(state["git_status_at_start"])
    for name in ("SUCCESS", "validation.json"):
        (report / name).unlink(missing_ok=True)
    return state


def finish_stage(report: Path, root: Path, state: dict, before: list, config_path: Path,
                 config: dict, products: list[Path], upstream: dict, versions: dict, command: list[str]) -> None:
    if sorted(file_records([root / row["path"] for row in before], root), key=lambda r: r["path"]) != sorted(before, key=lambda r: r["path"]):
        raise ValueError("Inputs changed during execution")
    shutil.copyfile(config_path, report / "config_snapshot.yaml")
    write_json(report / "input_checksums.json", before)
    outputs = [*products, report / "config_snapshot.yaml", report / "input_checksums.json"]
    write_json(report / "provenance.json", {**state, "status": "success", "parameters": config,
        "finished_at_utc": datetime.now(timezone.utc).isoformat(), "checksum_algorithm": "SHA-256",
        "inputs": before, "outputs": file_records(outputs, root), "upstream_validation": upstream,
        "software_versions": versions, "command": command})
    (report / "SUCCESS").write_text(state["git_commit"] + "\n", encoding="ascii")


def validate_artifacts(report: Path, root: Path, expected: list[Path], config: dict) -> dict:
    record = json.loads((report / "provenance.json").read_text())
    revision = (report / "SUCCESS").read_text().strip()
    if not re.fullmatch(r"[0-9a-f]{40}", revision) or record["status"] != "success" or record["git_commit"] != revision:
        raise ValueError("Invalid success or provenance record")
    for section in ("inputs", "outputs"):
        rows = record[section]
        paths = [repository_path(root, row["path"]) for row in rows]
        if len(set(paths)) != len(paths) or not rows:
            raise ValueError("Missing or duplicate checksum paths")
        if section == "outputs" and set(paths) != set(path.resolve() for path in expected):
            raise ValueError("Output membership mismatch")
        if sorted(file_records(paths, root), key=lambda r: r["path"]) != sorted(rows, key=lambda r: r["path"]):
            raise ValueError(f"{section} checksum mismatch")
    if json.loads((report / "input_checksums.json").read_text()) != record["inputs"]:
        raise ValueError("Input checksum list disagrees with provenance")
    if json.loads((report / "config_snapshot.yaml").read_text()) != config or record["parameters"] != config:
        raise ValueError("Configuration mismatch")
    return {"status": "valid", "source_git_commit": revision,
            "inputs_checked": len(record["inputs"]), "artifacts_checked": len(record["outputs"])}
