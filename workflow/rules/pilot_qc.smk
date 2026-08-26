"""Raw sequencing QC for the configured pilot vertical slice."""

import csv


PILOT_FASTQ_BY_READ = {}
PILOT_MANIFEST_BY_RUN = {}
with open("metadata/pilot_manifest.tsv", encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        PILOT_MANIFEST_BY_RUN[row["run_accession"]] = row
        for direction in ("1", "2"):
            read_id = f"{row['run_accession']}_{direction}"
            PILOT_FASTQ_BY_READ[read_id] = (
                f"data/raw/{row['study_id']}/{row['run_accession']}/"
                f"{read_id}.fastq.gz"
            )

PILOT_READ_IDS = sorted(PILOT_FASTQ_BY_READ)
PILOT_RUNS = sorted(PILOT_MANIFEST_BY_RUN)
PILOT_STUDY_IDS = {row["study_id"] for row in PILOT_MANIFEST_BY_RUN.values()}
if len(PILOT_STUDY_IDS) != 1:
    raise ValueError("Pilot sequencing workflow requires exactly one study")
PILOT_STUDY_ID = next(iter(PILOT_STUDY_IDS))
RAW_FASTQC_DIRECTORY = "results/qc/pilot/fastqc_raw"
RAW_MULTIQC_DIRECTORY = "results/qc/pilot/multiqc_raw"
CUTADAPT_DIRECTORY = "results/qc/pilot/cutadapt"
PILOT_INTERIM_PATTERN = (
    f"data/interim/{PILOT_STUDY_ID}/pilot/"
    "{run_accession}/{run_accession}_{direction}.fastq.gz"
)
PILOT_PROCESSING_CONFIG = config["amplicon_processing"]["PRJNA492720"]
PILOT_PRIMERS = PILOT_PROCESSING_CONFIG["primers"]
PILOT_CUTADAPT = PILOT_PROCESSING_CONFIG["cutadapt"]


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


rule pilot_primer_trimming:
    input:
        summary="results/qc/pilot/cutadapt_summary.tsv",
        trimmed=expand(
            PILOT_INTERIM_PATTERN,
            run_accession=PILOT_RUNS,
            direction=("1", "2"),
        ),


rule trim_pilot_primers:
    input:
        r1=lambda wildcards: PILOT_FASTQ_BY_READ[f"{wildcards.run_accession}_1"],
        r2=lambda wildcards: PILOT_FASTQ_BY_READ[f"{wildcards.run_accession}_2"],
    output:
        r1=PILOT_INTERIM_PATTERN.replace("{direction}", "1"),
        r2=PILOT_INTERIM_PATTERN.replace("{direction}", "2"),
        report=f"{CUTADAPT_DIRECTORY}/{{run_accession}}.cutadapt.json",
    params:
        output_directory=lambda wildcards: (
            f"data/interim/{PILOT_STUDY_ID}/pilot/{wildcards.run_accession}"
        ),
        report_directory=CUTADAPT_DIRECTORY,
        forward_adapter=(
            f"515F={PILOT_PRIMERS['forward']['gene_specific_sequence']};"
            f"min_overlap={PILOT_CUTADAPT['forward_minimum_overlap']}"
        ),
        reverse_adapter=(
            f"806R={PILOT_PRIMERS['reverse']['gene_specific_sequence']};"
            f"min_overlap={PILOT_CUTADAPT['reverse_minimum_overlap']}"
        ),
        error_rate=PILOT_CUTADAPT["error_rate"],
    threads: 2
    shell:
        "mkdir -p {params.output_directory:q} {params.report_directory:q} && "
        "cutadapt --cores {threads} --no-indels --error-rate {params.error_rate} "
        "--front {params.forward_adapter:q} -G {params.reverse_adapter:q} "
        "--json {output.report:q} --output {output.r1:q} "
        "--paired-output {output.r2:q} {input.r1:q} {input.r2:q}"


rule summarize_pilot_cutadapt:
    input:
        config="config/config.yaml",
        manifest="metadata/pilot_manifest.tsv",
        reports=expand(
            f"{CUTADAPT_DIRECTORY}/{{run_accession}}.cutadapt.json",
            run_accession=PILOT_RUNS,
        ),
        trimmed=expand(
            PILOT_INTERIM_PATTERN,
            run_accession=PILOT_RUNS,
            direction=("1", "2"),
        ),
    output:
        "results/qc/pilot/cutadapt_summary.tsv"
    params:
        report_directory=CUTADAPT_DIRECTORY,
    shell:
        "python scripts/qc/summarize_cutadapt.py --config {input.config:q} "
        "--manifest {input.manifest:q} "
        "--report-directory {params.report_directory:q} --output {output:q}"


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
