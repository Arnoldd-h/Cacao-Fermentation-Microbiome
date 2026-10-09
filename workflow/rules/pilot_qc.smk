"""Manifest-driven pilot QC; no data files are opened while parsing the workflow."""

from cacao_inventory.pilot_workflow import load_pilot_manifest, manifest_row, fastq_path, processing_configuration
from cacao_inventory.provenance import sha256_file


def pilot_manifest_path(wildcards=None):
    return checkpoints.build_pilot_manifest.get().output.manifest


def pilot_manifest_input(wildcards=None):
    return ancient(str(pilot_manifest_path(wildcards)))


def pilot_manifest_fingerprint(wildcards=None):
    return sha256_file(pilot_manifest_path(wildcards))


def pilot_manifest_rows(wildcards=None):
    return load_pilot_manifest(pilot_manifest_path(wildcards))


def pilot_study_id(wildcards=None):
    return pilot_manifest_rows(wildcards)[0]["study_id"]


def pilot_run_ids(wildcards=None):
    return sorted(row["run_accession"] for row in pilot_manifest_rows(wildcards))


def pilot_fastq_paths(wildcards=None, stage="trimmed"):
    return [fastq_path(row, direction, stage) for row in pilot_manifest_rows(wildcards) for direction in ("1", "2")]


def pilot_raw_fastq_paths(wildcards):
    return pilot_fastq_paths(wildcards, stage="raw")


def pilot_read_path(wildcards, stage="raw", direction=None):
    row = manifest_row(pilot_manifest_rows(wildcards), wildcards.run_accession, getattr(wildcards, "study_id", None))
    return fastq_path(row, direction or wildcards.direction, stage)


def pilot_processing(wildcards):
    processing = processing_configuration(config, pilot_manifest_rows(wildcards)[0])
    # Changing downstream DADA2 parameters must not repeat primer trimming.
    return {key: processing[key] for key in ("study_id", "primers", "primer_detection", "cutadapt")}


def pilot_fastqc_archives(wildcards, stage="raw"):
    return [f"results/qc/pilot/fastqc_{stage}/{run}_{direction}_fastqc.zip" for run in pilot_run_ids(wildcards) for direction in ("1", "2")]


QC_ENVIRONMENT = "environment/conda-linux-64.lock"
QC_COMMON_CODE = ["python/cacao_inventory/qc.py", "python/cacao_inventory/io.py", "python/cacao_inventory/schema.py", "python/cacao_inventory/config.py", "python/cacao_inventory/download.py", "python/cacao_inventory/provenance.py", "python/cacao_inventory/pilot_workflow.py", QC_ENVIRONMENT]


wildcard_constraints:
    direction="[12]",
    run_accession="(?:SRR|ERR|DRR)[0-9]+",
    study_id="[A-Za-z0-9][A-Za-z0-9_-]*",


rule download_pilot_read:
    input:
        # Raw data are immutable. The aggregate report revalidates bytes and MD5;
        # configuration changes cannot overwrite a previously downloaded file.
        manifest=pilot_manifest_input,
        config=ancient("config/config.yaml"),
        code=ancient(python_sources("scripts/qc/download_pilot_read.py", "download", "pilot_workflow", "config", "io", "schema")),
    output:
        protected("data/raw/{study_id}/{run_accession}/{run_accession}_{direction}.fastq.gz")
    shell:
        "python scripts/qc/download_pilot_read.py --manifest {input.manifest:q} --config {input.config:q} "
        "--study-id {wildcards.study_id:q} --run-accession {wildcards.run_accession:q} "
        "--direction {wildcards.direction:q} --output {output:q}"


rule download_pilot_fastq:
    input:
        manifest=pilot_manifest_input,
        fastq=pilot_raw_fastq_paths,
        config=ancient("config/config.yaml"),
        code=python_sources("scripts/metadata/download_pilot_fastq.py", "download", "config", "io", "schema"),
    output:
        "results/qc/pilot_download_validation.tsv"
    params:
        settings=config["pilot_download"],
        manifest_sha256=pilot_manifest_fingerprint,
    shell:
        "python scripts/metadata/download_pilot_fastq.py --manifest {input.manifest:q} --report {output:q}"


rule validate_pilot_fastq:
    input:
        report="results/qc/pilot_download_validation.tsv",
        fastq=pilot_raw_fastq_paths,
        code=python_sources("scripts/metadata/validate_pilot_fastq.py", "fastq", "io", "schema"),
    output:
        "results/qc/pilot_fastq_validation.tsv"
    shell:
        "python scripts/metadata/validate_pilot_fastq.py --download-report {input.report:q} --output {output:q}"


rule pilot_raw_qc:
    input:
        "results/qc/pilot/multiqc_raw/multiqc_report.html",
        "results/qc/pilot/raw_read_quality.tsv",
        "results/qc/pilot/raw_read_quality.provenance.json",
        manifest=pilot_manifest_input,


rule pilot_primer_detection:
    input:
        "results/qc/pilot/primer_detection.tsv",
        "results/qc/pilot/primer_detection.provenance.json",
        manifest=pilot_manifest_input,


rule detect_pilot_primers:
    input:
        config=ancient("config/config.yaml"),
        manifest=pilot_manifest_input,
        validated="results/qc/pilot_fastq_validation.tsv",
        fastq=pilot_raw_fastq_paths,
        code=["scripts/qc/detect_primers.py", *QC_COMMON_CODE],
    output:
        table="results/qc/pilot/primer_detection.tsv",
        provenance="results/qc/pilot/primer_detection.provenance.json",
    params:
        processing=pilot_processing,
        manifest_sha256=pilot_manifest_fingerprint,
    shell:
        "python scripts/qc/detect_primers.py --config {input.config:q} --manifest {input.manifest:q} --output {output.table:q}"


rule pilot_primer_trimming:
    input:
        "results/qc/pilot/cutadapt_summary.tsv",
        "results/qc/pilot/cutadapt_summary.provenance.json",
        pilot_fastq_paths,


rule trim_pilot_primers:
    input:
        manifest=pilot_manifest_input,
        config=ancient("config/config.yaml"),
        r1=lambda wildcards: pilot_read_path(wildcards, direction="1"),
        r2=lambda wildcards: pilot_read_path(wildcards, direction="2"),
        validated="results/qc/pilot_fastq_validation.tsv",
        detection="results/qc/pilot/primer_detection.tsv",
        code=["scripts/qc/trim_primers.py", *QC_COMMON_CODE],
    output:
        r1="data/interim/{study_id}/pilot/{run_accession}/{run_accession}_1.fastq.gz",
        r2="data/interim/{study_id}/pilot/{run_accession}/{run_accession}_2.fastq.gz",
        report="results/qc/pilot/cutadapt/{study_id}/{run_accession}.cutadapt.json",
    params:
        processing=pilot_processing,
        manifest_sha256=pilot_manifest_fingerprint,
    threads: 2
    shell:
        "python scripts/qc/trim_primers.py --config {input.config:q} --manifest {input.manifest:q} "
        "--run-accession {wildcards.run_accession:q} --study-id {wildcards.study_id:q} --threads {threads} "
        "--report {output.report:q}"


rule summarize_pilot_cutadapt:
    input:
        config=ancient("config/config.yaml"),
        manifest=pilot_manifest_input,
        reports=lambda wildcards: [f"results/qc/pilot/cutadapt/{pilot_study_id(wildcards)}/{run}.cutadapt.json" for run in pilot_run_ids(wildcards)],
        trimmed=pilot_fastq_paths,
        code=["scripts/qc/summarize_cutadapt.py", *QC_COMMON_CODE],
    output:
        table="results/qc/pilot/cutadapt_summary.tsv",
        provenance="results/qc/pilot/cutadapt_summary.provenance.json",
    params:
        report_directory=lambda wildcards: f"results/qc/pilot/cutadapt/{pilot_study_id(wildcards)}",
        processing=pilot_processing,
        manifest_sha256=pilot_manifest_fingerprint,
    shell:
        "python scripts/qc/summarize_cutadapt.py --config {input.config:q} --manifest {input.manifest:q} "
        "--report-directory {params.report_directory:q} --output {output.table:q}"


rule pilot_post_trim_qc:
    input:
        "results/qc/pilot/multiqc_trimmed/multiqc_report.html",
        "results/qc/pilot/trimmed_read_quality.tsv",
        "results/qc/pilot/read_quality_comparison.tsv",
        "results/qc/pilot/primer_detection_trimmed.tsv",
        "results/qc/pilot/trimmed_read_quality.provenance.json",
        "results/qc/pilot/primer_detection_trimmed.provenance.json",
        manifest=pilot_manifest_input,


rule fastqc_raw:
    input:
        fastq=lambda wildcards: pilot_read_path(wildcards),
        validated="results/qc/pilot_fastq_validation.tsv",
        environment=QC_ENVIRONMENT,
    output:
        html="results/qc/pilot/fastqc_raw/{run_accession}_{direction}_fastqc.html",
        archive="results/qc/pilot/fastqc_raw/{run_accession}_{direction}_fastqc.zip",
    threads: 2
    shell:
        "mkdir -p results/qc/pilot/fastqc_raw && fastqc --threads {threads} --outdir results/qc/pilot/fastqc_raw {input.fastq:q}"


rule fastqc_trimmed:
    input:
        fastq=lambda wildcards: pilot_read_path(wildcards, stage="trimmed"),
        environment=QC_ENVIRONMENT,
    output:
        html="results/qc/pilot/fastqc_trimmed/{run_accession}_{direction}_fastqc.html",
        archive="results/qc/pilot/fastqc_trimmed/{run_accession}_{direction}_fastqc.zip",
    threads: 2
    shell:
        "mkdir -p results/qc/pilot/fastqc_trimmed && fastqc --threads {threads} --outdir results/qc/pilot/fastqc_trimmed {input.fastq:q}"


rule multiqc_raw:
    input:
        archives=pilot_fastqc_archives,
        environment=QC_ENVIRONMENT,
    output:
        report="results/qc/pilot/multiqc_raw/multiqc_report.html",
        fastqc_table="results/qc/pilot/multiqc_raw/multiqc_report_data/multiqc_fastqc.txt",
    shell:
        "mkdir -p results/qc/pilot/multiqc_raw && multiqc --force --outdir results/qc/pilot/multiqc_raw "
        "--filename multiqc_report.html {input.archives:q}"


rule multiqc_trimmed:
    input:
        archives=lambda wildcards: pilot_fastqc_archives(wildcards, stage="trimmed"),
        environment=QC_ENVIRONMENT,
    output:
        report="results/qc/pilot/multiqc_trimmed/multiqc_report.html",
        fastqc_table="results/qc/pilot/multiqc_trimmed/multiqc_report_data/multiqc_fastqc.txt",
    shell:
        "mkdir -p results/qc/pilot/multiqc_trimmed && multiqc --force --outdir results/qc/pilot/multiqc_trimmed "
        "--filename multiqc_report.html {input.archives:q}"


rule summarize_raw_fastqc:
    input:
        manifest=pilot_manifest_input,
        multiqc="results/qc/pilot/multiqc_raw/multiqc_report_data/multiqc_fastqc.txt",
        validation="results/qc/pilot_fastq_validation.tsv",
        config=ancient("config/config.yaml"),
        code=["scripts/qc/summarize_fastqc.py", *QC_COMMON_CODE],
    output:
        table="results/qc/pilot/raw_read_quality.tsv",
        provenance="results/qc/pilot/raw_read_quality.provenance.json",
    params:
        manifest_sha256=pilot_manifest_fingerprint,
    shell:
        "python scripts/qc/summarize_fastqc.py --multiqc-fastqc {input.multiqc:q} "
        "--fastq-validation {input.validation:q} --output {output.table:q}"


rule summarize_trimmed_fastqc:
    input:
        manifest=pilot_manifest_input,
        multiqc="results/qc/pilot/multiqc_trimmed/multiqc_report_data/multiqc_fastqc.txt",
        cutadapt="results/qc/pilot/cutadapt_summary.tsv",
        raw="results/qc/pilot/raw_read_quality.tsv",
        config=ancient("config/config.yaml"),
        code=["scripts/qc/summarize_trimmed_fastqc.py", *QC_COMMON_CODE],
    output:
        metrics="results/qc/pilot/trimmed_read_quality.tsv",
        comparison="results/qc/pilot/read_quality_comparison.tsv",
        provenance="results/qc/pilot/trimmed_read_quality.provenance.json",
    params:
        manifest_sha256=pilot_manifest_fingerprint,
    shell:
        "python scripts/qc/summarize_trimmed_fastqc.py --multiqc-fastqc {input.multiqc:q} "
        "--cutadapt-summary {input.cutadapt:q} --raw-quality {input.raw:q} "
        "--output {output.metrics:q} --comparison-output {output.comparison:q}"


rule detect_trimmed_primers:
    input:
        config=ancient("config/config.yaml"),
        manifest=pilot_manifest_input,
        fastq=pilot_fastq_paths,
        code=["scripts/qc/detect_primers.py", *QC_COMMON_CODE],
    output:
        table="results/qc/pilot/primer_detection_trimmed.tsv",
        provenance="results/qc/pilot/primer_detection_trimmed.provenance.json",
    params:
        processing=pilot_processing,
        manifest_sha256=pilot_manifest_fingerprint,
    shell:
        "python scripts/qc/detect_primers.py --input-stage trimmed --config {input.config:q} "
        "--manifest {input.manifest:q} --output {output.table:q}"
