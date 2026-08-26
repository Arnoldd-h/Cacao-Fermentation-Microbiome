"""Raw sequencing QC for the configured pilot vertical slice."""

import csv


PILOT_FASTQ_BY_READ = {}
with open("metadata/pilot_manifest.tsv", encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        for direction in ("1", "2"):
            read_id = f"{row['run_accession']}_{direction}"
            PILOT_FASTQ_BY_READ[read_id] = (
                f"data/raw/{row['study_id']}/{row['run_accession']}/"
                f"{read_id}.fastq.gz"
            )

PILOT_READ_IDS = sorted(PILOT_FASTQ_BY_READ)
RAW_FASTQC_DIRECTORY = "results/qc/pilot/fastqc_raw"
RAW_MULTIQC_DIRECTORY = "results/qc/pilot/multiqc_raw"


rule pilot_raw_qc:
    input:
        report=f"{RAW_MULTIQC_DIRECTORY}/multiqc_report.html",
        metrics="results/qc/pilot/raw_read_quality.tsv",


rule pilot_primer_detection:
    input:
        "results/qc/pilot/primer_detection.tsv",


rule detect_pilot_primers:
    input:
        config="config/config.yaml",
        manifest="metadata/pilot_manifest.tsv",
        fastq=list(PILOT_FASTQ_BY_READ.values()),
    output:
        "results/qc/pilot/primer_detection.tsv"
    shell:
        "python scripts/qc/detect_primers.py "
        "--config {input.config:q} --manifest {input.manifest:q} --output {output:q}"


rule fastqc_raw:
    input:
        lambda wildcards: PILOT_FASTQ_BY_READ[wildcards.read_id]
    output:
        html=f"{RAW_FASTQC_DIRECTORY}/{{read_id}}_fastqc.html",
        archive=f"{RAW_FASTQC_DIRECTORY}/{{read_id}}_fastqc.zip",
    params:
        output_directory=RAW_FASTQC_DIRECTORY,
    threads: 2
    shell:
        "mkdir -p {params.output_directory:q} && "
        "fastqc --threads {threads} --outdir {params.output_directory:q} {input:q}"


rule multiqc_raw:
    input:
        expand(f"{RAW_FASTQC_DIRECTORY}/{{read_id}}_fastqc.zip", read_id=PILOT_READ_IDS)
    output:
        report=f"{RAW_MULTIQC_DIRECTORY}/multiqc_report.html",
        fastqc_table=(
            f"{RAW_MULTIQC_DIRECTORY}/multiqc_report_data/multiqc_fastqc.txt"
        ),
    params:
        output_directory=RAW_MULTIQC_DIRECTORY,
        input_directory=RAW_FASTQC_DIRECTORY,
    shell:
        "mkdir -p {params.output_directory:q} && "
        "multiqc --force --outdir {params.output_directory:q} "
        "--filename multiqc_report.html {params.input_directory:q}"


rule summarize_raw_fastqc:
    input:
        multiqc=f"{RAW_MULTIQC_DIRECTORY}/multiqc_report_data/multiqc_fastqc.txt",
        validation="results/qc/pilot_fastq_validation.tsv",
    output:
        "results/qc/pilot/raw_read_quality.tsv"
    shell:
        "python scripts/qc/summarize_fastqc.py "
        "--multiqc-fastqc {input.multiqc:q} "
        "--fastq-validation {input.validation:q} --output {output:q}"
