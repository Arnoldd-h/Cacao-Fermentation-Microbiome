"""Pinned taxonomy of the independently validated DADA2 pilot."""

from cacao_inventory.config import load_json_yaml

taxonomy_config = load_json_yaml("config/taxonomy.yaml")
taxonomy_reference = taxonomy_config["reference"]
taxonomy_products = ["taxonomy.tsv", "taxonomy_unfiltered.tsv", "taxonomy_bootstraps.tsv",
    "taxonomy_sensitivity.tsv", "taxonomy_screening.tsv", "assignment_coverage.tsv", "sample_coverage.tsv",
    "software_versions.tsv", "session_info.txt", "warnings.txt", "config_snapshot.yaml", "input_checksums.json",
    "assignment_coverage.pdf", "assignment_coverage.svg", "assignment_coverage.png", "provenance.json", "SUCCESS"]


rule pilot_taxonomy:
    input:
        "results/taxonomy/pilot/validation.json",
        manifest=pilot_manifest_input,


rule download_taxonomy_database:
    input:
        config=ancient("config/taxonomy.yaml"),
        code=ancient(python_sources("scripts/taxonomy/download_reference.py", "config", "download", "provenance", "taxonomy_reference")),
    output:
        protected(taxonomy_reference["path"]),
    shell:
        "python scripts/taxonomy/download_reference.py --config {input.config:q} --database-only"


rule validate_taxonomy_reference:
    input:
        reference=taxonomy_reference["path"],
        config="config/taxonomy.yaml",
        code=python_sources("scripts/taxonomy/download_reference.py", "config", "download", "provenance", "taxonomy_reference"),
    output:
        provenance=taxonomy_reference["provenance_path"],
    params:
        config_sha256=sha256_file("config/taxonomy.yaml"),
    shell:
        "python scripts/taxonomy/download_reference.py --config {input.config:q} --validate-only"


rule run_pilot_taxonomy:
    input:
        manifest=pilot_manifest_input,
        dada2="results/dada2/pilot/validation.json",
        sequences="results/dada2/pilot/asv_sequences.tsv",
        counts="results/dada2/pilot/asv_counts.tsv",
        metadata="results/dada2/pilot/sample_metadata.tsv",
        reference=taxonomy_reference["path"],
        reference_provenance=taxonomy_reference["provenance_path"],
        config="config/taxonomy.yaml",
        project_config=ancient("config/config.yaml"),
        code=python_sources("scripts/taxonomy/run_taxonomy.py", "config", "download", "provenance", "taxonomy_reference") +
            ["scripts/taxonomy/assign_taxonomy.R", "scripts/taxonomy/helpers.R", "scripts/dada2/helpers.R", "scripts/dada2/validate_outputs.py",
             "workflow/Snakefile", "workflow/rules/pilot_taxonomy.smk"],
        environment="environment/conda-linux-64.lock",
    output:
        expand("results/taxonomy/pilot/{product}", product=taxonomy_products),
    params:
        config_sha256=sha256_file("config/taxonomy.yaml"),
        project_config_sha256=sha256_file("config/config.yaml"),
        reference_md5=taxonomy_reference["expected_md5"],
    threads: taxonomy_config["classification"]["threads"]
    log:
        "results/intermediate/taxonomy/pilot/execution.log"
    shell:
        "mkdir -p results/intermediate/taxonomy/pilot && "
        "python scripts/taxonomy/run_taxonomy.py --config {input.config:q} "
        "--threads {threads} > {log:q} 2>&1"


rule validate_pilot_taxonomy:
    input:
        manifest=pilot_manifest_input,
        products=rules.run_pilot_taxonomy.output,
        code=python_sources("scripts/taxonomy/validate_taxonomy.py", "config", "io", "provenance", "taxonomy_reference", "taxonomy_validation"),
    output:
        "results/taxonomy/pilot/validation.json"
    shell:
        "python scripts/taxonomy/validate_taxonomy.py --output {output:q}"
