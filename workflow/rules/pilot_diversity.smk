"""Auditable bacterial separation and descriptive within-study pilot diversity."""

from cacao_inventory.config import load_json_yaml

diversity_config = load_json_yaml("config/diversity.yaml")
diversity_paths = diversity_config["paths"]
bacterial_products = ["asv_counts.tsv", "asv_sequences.tsv", "taxonomy.tsv", "sample_metadata.tsv"]
filtering_products = ["asv_filter_log.tsv", "sample_retention.tsv", "summary.json", "config_snapshot.yaml",
                      "input_checksums.json", "provenance.json", "SUCCESS"]
diversity_products = ["alpha_diversity.tsv", "clr_coordinates.tsv", "aitchison_distances.tsv", "bray_curtis_distances.tsv",
    "pca_scores.tsv", "pca_variance.tsv", "software_versions.tsv", "session_info.txt", "warnings.txt",
    "config_snapshot.yaml", "input_checksums.json", "provenance.json", "SUCCESS",
    *[name + "." + extension for name in ("alpha_diversity", "aitchison_pca") for extension in ("pdf", "svg", "png")]]


rule pilot_bacterial_table:
    input:
        diversity_paths["filtering_dir"] + "/validation.json",
        manifest=pilot_manifest_input,


rule pilot_diversity:
    input:
        diversity_paths["diversity_dir"] + "/validation.json",
        manifest=pilot_manifest_input,


rule prepare_pilot_bacterial_table:
    input:
        manifest=pilot_manifest_input,
        dada2=diversity_paths["dada2_dir"] + "/validation.json",
        taxonomy=diversity_paths["taxonomy_dir"] + "/validation.json",
        config="config/diversity.yaml",
        sources=[diversity_paths["dada2_dir"] + "/" + name for name in ("asv_counts.tsv", "asv_sequences.tsv", "sample_metadata.tsv")] +
                [diversity_paths["taxonomy_dir"] + "/" + name for name in ("taxonomy.tsv", "config_snapshot.yaml", "provenance.json")],
        code=python_sources("scripts/filtering/prepare_bacterial_table.py", "analysis_artifacts", "bacterial_filter", "config", "io", "provenance", "taxonomy_reference", "taxonomy_validation") +
             ["scripts/dada2/validate_outputs.py", "scripts/taxonomy/validate_taxonomy.py", "workflow/rules/pilot_diversity.smk", "environment/conda-linux-64.lock"],
    output:
        [diversity_paths["bacterial_dir"] + "/" + name for name in bacterial_products] +
        [diversity_paths["filtering_dir"] + "/" + name for name in filtering_products],
    params:
        config_sha256=sha256_file("config/diversity.yaml"),
    shell:
        "python scripts/filtering/prepare_bacterial_table.py --config {input.config:q}"


rule validate_pilot_bacterial_table:
    input:
        manifest=pilot_manifest_input,
        products=rules.prepare_pilot_bacterial_table.output,
        config="config/diversity.yaml",
        code=python_sources("scripts/filtering/validate_bacterial_table.py", "analysis_artifacts", "bacterial_filter", "config", "io", "provenance", "taxonomy_reference"),
    output:
        diversity_paths["filtering_dir"] + "/validation.json",
    shell:
        "python scripts/filtering/validate_bacterial_table.py --config {input.config:q} --output {output:q}"


rule describe_pilot_diversity:
    input:
        manifest=pilot_manifest_input,
        filtering=diversity_paths["filtering_dir"] + "/validation.json",
        filtering_provenance=diversity_paths["filtering_dir"] + "/provenance.json",
        bacterial=[diversity_paths["bacterial_dir"] + "/" + name for name in bacterial_products],
        config="config/diversity.yaml",
        project_config=ancient("config/config.yaml"),
        code=python_sources("scripts/diversity/run_diversity.py", "analysis_artifacts", "bacterial_filter", "config", "io", "provenance", "taxonomy_reference") +
             ["scripts/diversity/describe_pilot.R", "scripts/dada2/helpers.R", "scripts/filtering/validate_bacterial_table.py",
              "workflow/rules/pilot_diversity.smk", "environment/conda-linux-64.lock"],
    output:
        [diversity_paths["diversity_dir"] + "/" + name for name in diversity_products],
    params:
        config_sha256=sha256_file("config/diversity.yaml"),
        project_sha256=sha256_file("config/config.yaml"),
    shell:
        "python scripts/diversity/run_diversity.py --config {input.config:q}"


rule validate_pilot_diversity:
    input:
        manifest=pilot_manifest_input,
        products=rules.describe_pilot_diversity.output,
        config="config/diversity.yaml",
        code=python_sources("scripts/diversity/validate_diversity.py", "analysis_artifacts", "diversity_validation", "config", "io", "provenance", "taxonomy_reference") +
             ["scripts/diversity/run_diversity.py"],
    output:
        diversity_paths["diversity_dir"] + "/validation.json",
    shell:
        "python scripts/diversity/validate_diversity.py --config {input.config:q} --output {output:q}"
