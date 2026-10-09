"""Independent formula and geometry checks of R descriptive diversity outputs."""

import math
from pathlib import Path

from .analysis_artifacts import read_counts, read_table, unique_rows


def close(value: str, expected: float, tolerance: float, label: str) -> None:
    observed = float(value)
    if not math.isfinite(observed) or not math.isclose(observed, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise ValueError(f"Invalid {label}: {observed} != {expected}")


def validate_diversity_tables(directory: Path, bacterial: Path, config: dict) -> dict:
    rows, features = read_counts(bacterial / "asv_counts.tsv")
    counts = unique_rows(rows, ("study_id", "sample_id"))
    study = rows[0]["study_id"]
    samples = [row["sample_id"] for row in rows]
    settings = config["analysis"]
    tolerance = settings["validation_tolerance"]
    metadata = unique_rows(read_table(bacterial / "sample_metadata.tsv", ["study_id", "sample_id"]), ("study_id", "sample_id"))
    alpha = unique_rows(read_table(directory / "alpha_diversity.tsv", ["study_id", "sample_id", "library_reads", *settings["alpha_metrics"]]), ("study_id", "sample_id"))
    if set(metadata) != set(counts) or set(alpha) != set(counts):
        raise ValueError("Diversity sample membership mismatch")
    proportions = {}
    for key, row in counts.items():
        integers = [int(row[feature]) for feature in features]
        total = sum(integers)
        if not total:
            raise ValueError("Empty sample")
        p = [value / total for value in integers]
        proportions[key[1]] = p
        entropy = -math.fsum(value * math.log(value) for value in p if value)
        concentration = math.fsum(value**2 for value in p)
        expected = {"library_reads": total, "observed_asvs": sum(value > 0 for value in integers),
                    "shannon": entropy, "gini_simpson": 1 - concentration, "inverse_simpson": 1 / concentration,
                    "shannon_effective": math.exp(entropy)}
        for field, value in expected.items():
            close(alpha[key][field], value, tolerance, field)
        for field in ("run_accession", "fermentation_batch", "fermentation_hours", "relative_time", "fermentation_stage", "sampling_stratum"):
            if alpha[key].get(field) != metadata[key][field]:
                raise ValueError("Diversity metadata changed")
    constants = [settings["primary_pseudocount"], *settings["sensitivity_pseudocounts"]]
    # Numeric normalization accepts R's '1' for a configured JSON 1.0.
    def keyed_numeric(table, fields):
        for row in table:
            if "pseudocount" in row:
                row["pseudocount"] = str(float(row["pseudocount"]))
        return unique_rows(table, fields)
    clr = keyed_numeric(read_table(directory / "clr_coordinates.tsv", ["study_id", "pseudocount", "sample_id", "asv_id", "clr"]),
                        ("study_id", "pseudocount", "sample_id", "asv_id"))
    expected_clr = {(study, str(float(c)), sample, feature) for c in constants for sample in samples for feature in features}
    if set(clr) != expected_clr:
        raise ValueError("CLR ASV/sample/pseudocount membership mismatch")
    computed = {}
    for constant in constants:
        label = str(float(constant))
        for sample in samples:
            logs = [math.log(int(counts[study, sample][feature]) + constant) for feature in features]
            mean = math.fsum(logs) / len(features)
            vector = [value - mean for value in logs]
            computed[label, sample] = vector
            for feature, value in zip(features, vector):
                close(clr[study, label, sample, feature]["clr"], value, tolerance, "CLR")
    aitchison = keyed_numeric(read_table(directory / "aitchison_distances.tsv", ["study_id", "pseudocount", "sample_id_1", "sample_id_2", "distance"]),
                              ("study_id", "pseudocount", "sample_id_1", "sample_id_2"))
    bray = unique_rows(read_table(directory / "bray_curtis_distances.tsv", ["study_id", "sample_id_1", "sample_id_2", "distance"]),
                       ("study_id", "sample_id_1", "sample_id_2"))
    if set(aitchison) != {(study, str(float(c)), a, b) for c in constants for a in samples for b in samples} or set(bray) != {(study, a, b) for a in samples for b in samples}:
        raise ValueError("Distance matrix membership mismatch")
    for a in samples:
        for b in samples:
            close(bray[study, a, b]["distance"], 0.5 * math.fsum(abs(x - y) for x, y in zip(proportions[a], proportions[b])), tolerance, "Bray-Curtis")
            for constant in constants:
                label = str(float(constant))
                distance = math.sqrt(math.fsum((x - y)**2 for x, y in zip(computed[label, a], computed[label, b])))
                close(aitchison[study, label, a, b]["distance"], distance, tolerance, "Aitchison")
    axes = list(range(1, min(len(samples) - 1, len(features)) + 1))
    scores = keyed_numeric(read_table(directory / "pca_scores.tsv", ["study_id", "pseudocount", "sample_id", "axis", "score"]),
                           ("study_id", "pseudocount", "sample_id", "axis"))
    variances = keyed_numeric(read_table(directory / "pca_variance.tsv", ["study_id", "pseudocount", "axis", "variance", "explained_fraction"]),
                              ("study_id", "pseudocount", "axis"))
    if set(scores) != {(study, str(float(c)), sample, str(axis)) for c in constants for sample in samples for axis in axes} or set(variances) != {(study, str(float(c)), str(axis)) for c in constants for axis in axes}:
        raise ValueError("PCA component membership mismatch")
    for constant in constants:
        label = str(float(constant))
        vectors = {sample: [float(scores[study, label, sample, str(axis)]["score"]) for axis in axes] for sample in samples}
        if any(not math.isfinite(value) for vector in vectors.values() for value in vector):
            raise ValueError("Nonfinite PCA scores")
        means = [math.fsum(computed[label, sample][j] for sample in samples) / len(samples) for j in range(len(features))]
        total_variance = math.fsum((computed[label, sample][j] - means[j])**2 for sample in samples for j in range(len(features))) / (len(samples) - 1)
        if total_variance <= 0:
            raise ValueError("PCA has no variation")
        previous = math.inf
        for j, axis in enumerate(axes):
            values = [vectors[sample][j] for sample in samples]
            close(str(math.fsum(values)), 0, tolerance * max(1, math.sqrt(total_variance)), "PCA centering")
            variance = math.fsum(value**2 for value in values) / (len(samples) - 1)
            row = variances[study, label, str(axis)]
            close(row["variance"], variance, tolerance, "PCA variance")
            close(row["explained_fraction"], variance / total_variance, tolerance, "PCA explained fraction")
            if variance > previous + tolerance * max(1, previous if math.isfinite(previous) else 1):
                raise ValueError("PCA variance order invalid")
            previous = variance
            for k in range(j):
                close(str(math.fsum(vectors[sample][j] * vectors[sample][k] for sample in samples)), 0,
                      tolerance * max(1, total_variance * len(samples)), "PCA orthogonality")
        close(str(math.fsum(float(variances[study, label, str(axis)]["explained_fraction"]) for axis in axes)), 1, tolerance, "PCA total variance")
        for a in samples:
            for b in samples:
                distance = math.sqrt(math.fsum((x - y)**2 for x, y in zip(vectors[a], vectors[b])))
                close(aitchison[study, label, a, b]["distance"], distance, tolerance, "PCA distance preservation")
    depths = [int(alpha[key]["library_reads"]) for key in counts]
    richness = [int(alpha[key]["observed_asvs"]) for key in counts]
    return {"study_id": study, "samples": len(samples), "asvs": len(features), "retained_reads": sum(depths),
            "library_reads_min": min(depths), "library_reads_max": max(depths),
            "observed_asvs_min": min(richness), "observed_asvs_max": max(richness),
            "independent_batches": len({row["fermentation_batch"] for row in metadata.values() if row["fermentation_batch"]}),
            "pseudocounts": constants, "rarefaction": False, "inference": False}
