"""DADA2 technical pilot, with explicit inputs and checked portable outputs."""


def pilot_dada2_settings(wildcards):
    row = pilot_manifest_rows(wildcards)[0]
    processing = processing_configuration(config, row)
    if "dada2" not in processing:
        raise ValueError(f"DADA2 parameters have not been registered for {row['bioproject']}")
    return processing["dada2"]


rule pilot_dada2:
    input:
        "results/dada2/pilot/validation.json",
        manifest=pilot_manifest_input,


rule run_pilot_dada2:
    input:
        manifest=pilot_manifest_input,
        config=ancient("config/config.yaml"),
        fastq=pilot_fastq_paths,
        raw_quality="results/qc/pilot/raw_read_quality.tsv",
        trimmed_quality="results/qc/pilot/trimmed_read_quality.tsv",
        primer_detection="results/qc/pilot/primer_detection_trimmed.tsv",
        qc_provenance=["results/qc/pilot/trimmed_read_quality.provenance.json", "results/qc/pilot/primer_detection_trimmed.provenance.json"],
        code=["scripts/dada2/run_pilot.R", "scripts/dada2/helpers.R"],
        environment="environment/conda-linux-64.lock",
    output:
        success="results/dada2/pilot/SUCCESS",
        tracking="results/dada2/pilot/read_tracking.tsv",
        counts="results/dada2/pilot/asv_counts.tsv",
        sequences="results/dada2/pilot/asv_sequences.tsv",
        fasta="results/dada2/pilot/asv_sequences.fasta",
        summary="results/dada2/pilot/summary.tsv",
        lengths="results/dada2/pilot/sequence_length_distribution.tsv",
        merge="results/dada2/pilot/merge_diagnostics.tsv",
        errors="results/dada2/pilot/error_learning.tsv",
        provenance="results/dada2/pilot/run_provenance.yaml",
        input_checksums="results/dada2/pilot/input_checksums.tsv",
        output_checksums="results/dada2/pilot/output_checksums.tsv",
        snapshot="results/dada2/pilot/config_snapshot.yaml",
        metadata="results/dada2/pilot/sample_metadata.tsv",
        versions="results/dada2/pilot/software_versions.tsv",
        session="results/dada2/pilot/session_info.txt",
        warnings="results/dada2/pilot/warnings.txt",
        figures=expand("results/dada2/pilot/error_model_{direction}.{extension}", direction=["F", "R"], extension=["pdf", "svg", "png"]),
    params:
        settings=pilot_dada2_settings,
        seed=config["project"]["default_random_seed"],
        manifest_sha256=pilot_manifest_fingerprint,
        config_sha256=sha256_file("config/config.yaml"),
    threads: lambda wildcards: int(pilot_dada2_settings(wildcards)["threads"])
    log:
        "results/intermediate/dada2/pilot/execution.log"
    shell:
        "mkdir -p results/intermediate/dada2/pilot && "
        "Rscript scripts/dada2/run_pilot.R --config {input.config:q} --manifest {input.manifest:q} "
        "--output-dir results/dada2/pilot --intermediate-dir results/intermediate/dada2/pilot "
        "--raw-quality {input.raw_quality:q} --trimmed-quality {input.trimmed_quality:q} "
        "--threads {threads} > {log:q} 2>&1"


rule validate_pilot_dada2:
    input:
        manifest=pilot_manifest_input,
        products=rules.run_pilot_dada2.output,
        code="scripts/dada2/validate_outputs.py",
    output:
        "results/dada2/pilot/validation.json"
    shell:
        "python scripts/dada2/validate_outputs.py --run-dir results/dada2/pilot --check-inputs-root . --output {output:q}"
