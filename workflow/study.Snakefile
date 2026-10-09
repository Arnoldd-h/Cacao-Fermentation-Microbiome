"""Complete single-study processing with unchanged pilot scientific methods."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd() / "python"))
from cacao_inventory.config import load_json_yaml
from cacao_inventory.pilot_workflow import load_pilot_manifest, fastq_path, manifest_row
from cacao_inventory.provenance import sha256_file
from cacao_inventory.study_qc import trimmed_path
from cacao_inventory.study_scope import validate_scope

configfile: "config/full_study.yaml"
validate_scope(config, Path.cwd())
scope = config
project = load_json_yaml(scope["project_config"])
diversity = load_json_yaml(scope["diversity_config"])["paths"]
qc = scope["qc_dir"]
interim = scope["interim_dir"]
dada = scope["dada2_dir"]
tax = scope["taxonomy_dir"]
screen_dir = "results/intermediate/qc/" + scope["study_id"]
new_code = ["scripts/study/qc_stage.py", "scripts/qc/detect_primers.py", "environment/conda-linux-64.lock",
            *["python/cacao_inventory/" + name + ".py" for name in ("study_qc", "study_scope", "qc", "config", "io", "schema", "provenance", "analysis_artifacts", "pilot_workflow")]]
processing = project["amplicon_processing"][scope["bioproject"]]
reference = load_json_yaml(scope["taxonomy_config"])["reference"]
analysis_code = [str(path) for path in Path("python/cacao_inventory").glob("*.py")]

def manifest_input(wildcards=None):
    return ancient(str(checkpoints.study_manifest.get().output.manifest))

def rows(wildcards=None):
    return load_pilot_manifest(str(checkpoints.study_manifest.get().output.manifest))

def manifest_hash(wildcards=None):
    return sha256_file(str(checkpoints.study_manifest.get().output.manifest))

def reads(wildcards=None, stage="raw"):
    return [fastq_path(row, d) if stage == "raw" else Path(trimmed_path(scope, row, d)).as_posix() for row in rows() for d in ("1", "2")]

def read_path(wildcards, stage="raw", direction=None):
    row = manifest_row(rows(), wildcards.run_accession, scope["study_id"])
    d = direction or wildcards.direction
    return fastq_path(row, d) if stage == "raw" else Path(trimmed_path(scope, row, d)).as_posix()

def archives(wildcards, stage="raw"):
    return [qc + "/fastqc_" + stage + "/" + row["run_accession"] + "_" + d + "_fastqc.zip" for row in rows() for d in ("1", "2")]

rule all:
    input:
        diversity["diversity_dir"] + "/validation.json",
        manifest=manifest_input,

rule study_qc:
    input: qc + "/validation.json", manifest=manifest_input

rule study_dada2:
    input: dada + "/validation.json", manifest=manifest_input

checkpoint study_manifest:
    input:
        config="config/full_study.yaml",
        sources=[scope["project_config"], scope["datasets_config"], *[scope["inventory_dir"] + "/" + n + ".tsv" for n in ("studies", "runs", "samples", "exclusion_log")]],
        code=["scripts/study/build_manifest.py", *analysis_code],
    output:
        manifest=scope["manifest"],
        products=[scope["plan_dir"] + "/" + n for n in ("source_audit.json", "resource_plan.json", "fastq_stream_log.tsv", "provenance.json", "input_checksums.json", "config_snapshot.yaml", "SUCCESS")],
    params: config_sha256=sha256_file("config/full_study.yaml")
    shell: "python scripts/study/build_manifest.py --config {input.config:q} --verify-sources"

wildcard_constraints:
    direction="[12]", run_accession="(?:SRR|ERR|DRR)[0-9]+", stage="raw|trimmed"

rule download_read:
    input: manifest=manifest_input, config=ancient(scope["project_config"])
    output: protected("data/raw/" + scope["study_id"] + "/{run_accession}/{run_accession}_{direction}.fastq.gz")
    params: study_id=scope["study_id"]
    shell:
        "python scripts/qc/download_pilot_read.py --manifest {input.manifest:q} --config {input.config:q} "
        "--study-id {params.study_id:q} --run-accession {wildcards.run_accession:q} --direction {wildcards.direction:q} --output {output:q}"

rule download_validation:
    input: manifest=manifest_input, reads=reads, config=ancient(scope["project_config"])
    output: qc + "/download_validation.tsv"
    params: manifest_sha256=manifest_hash
    shell: "python scripts/metadata/download_pilot_fastq.py --config {input.config:q} --manifest {input.manifest:q} --report {output:q}"

rule validate_fastq:
    input: report=rules.download_validation.output, reads=reads
    output: qc + "/fastq_validation.tsv"
    shell: "python scripts/metadata/validate_pilot_fastq.py --download-report {input.report:q} --output {output:q}"

rule fastqc:
    input:
        fastq=lambda wc: read_path(wc, wc.stage),
        validated=lambda wc: qc + "/fastq_validation.tsv",
        manifest=manifest_input,
    output:
        html=qc + "/fastqc_{stage}/{run_accession}_{direction}_fastqc.html",
        archive=qc + "/fastqc_{stage}/{run_accession}_{direction}_fastqc.zip",
    params: directory=lambda wc: qc + "/fastqc_" + wc.stage
    resources: mem_mb=scope["resources"]["fastqc_memory_mb"]
    shell: "mkdir -p {params.directory:q} && fastqc --threads 1 --outdir {params.directory:q} {input.fastq:q}"

rule multiqc:
    input: archives=lambda wc: archives(wc, wc.stage), manifest=manifest_input
    output:
        html=qc + "/multiqc_{stage}/multiqc_report.html",
        table=qc + "/multiqc_{stage}/multiqc_report_data/multiqc_fastqc.txt",
    params: directory=lambda wc: qc + "/multiqc_" + wc.stage
    shell: "mkdir -p {params.directory:q} && multiqc --force --outdir {params.directory:q} --filename multiqc_report.html {input.archives:q}"

rule raw_summary:
    input: table=qc + "/multiqc_raw/multiqc_report_data/multiqc_fastqc.txt", validation=qc + "/fastq_validation.tsv", manifest=manifest_input
    output: table=qc + "/raw_read_quality.tsv", provenance=qc + "/raw_read_quality.provenance.json"
    shell: "python scripts/qc/summarize_fastqc.py --multiqc-fastqc {input.table:q} --fastq-validation {input.validation:q} --output {output.table:q}"

rule primer_screen:
    input:
        manifest=manifest_input, config="config/full_study.yaml", project=ancient(scope["project_config"]), code=new_code,
        r1=lambda wc: read_path(wc, wc.stage, "1"), r2=lambda wc: read_path(wc, wc.stage, "2"),
        validation=qc + "/fastq_validation.tsv",
    output: screen_dir + "/primers_{stage}/{run_accession}.tsv"
    params: method=processing["primer_detection"], primers=processing["primers"], manifest_sha256=manifest_hash
    shell: "python scripts/study/qc_stage.py screen --config {input.config:q} --run-accession {wildcards.run_accession:q} --stage {wildcards.stage:q} --output {output:q}"

rule screen_summary:
    input:
        tables=lambda wc: [screen_dir + "/primers_" + wc.stage + "/" + row["run_accession"] + ".tsv" for row in rows()],
        manifest=manifest_input, config="config/full_study.yaml", code=new_code,
    output: qc + "/primer_detection_{stage}.tsv"
    shell: "python scripts/study/qc_stage.py summarize_screen --config {input.config:q} --stage {wildcards.stage:q}"

rule trim:
    input:
        manifest=manifest_input, config="config/full_study.yaml", project=ancient(scope["project_config"]), code=new_code,
        r1=lambda wc: read_path(wc, direction="1"), r2=lambda wc: read_path(wc, direction="2"),
        detection=qc + "/primer_detection_raw.tsv",
    output:
        r1=interim + "/{run_accession}/{run_accession}_1.fastq.gz",
        r2=interim + "/{run_accession}/{run_accession}_2.fastq.gz",
        report=qc + "/cutadapt/{run_accession}.cutadapt.json",
    params: method=processing["cutadapt"], primers=processing["primers"], manifest_sha256=manifest_hash
    shell: "python scripts/study/qc_stage.py trim --config {input.config:q} --run-accession {wildcards.run_accession:q}"

rule cutadapt_summary:
    input:
        reports=lambda wc: [qc + "/cutadapt/" + row["run_accession"] + ".cutadapt.json" for row in rows()],
        reads=lambda wc: reads(wc, "trimmed"), manifest=manifest_input, config="config/full_study.yaml", code=new_code,
    output: qc + "/cutadapt_summary.tsv"
    shell: "python scripts/study/qc_stage.py summarize_cutadapt --config {input.config:q}"

rule trimmed_summary:
    input:
        table=qc + "/multiqc_trimmed/multiqc_report_data/multiqc_fastqc.txt", cutadapt=qc + "/cutadapt_summary.tsv",
        raw=qc + "/raw_read_quality.tsv", manifest=manifest_input,
    output:
        table=qc + "/trimmed_read_quality.tsv", comparison=qc + "/read_quality_comparison.tsv", provenance=qc + "/trimmed_read_quality.provenance.json",
    shell:
        "python scripts/qc/summarize_trimmed_fastqc.py --multiqc-fastqc {input.table:q} --cutadapt-summary {input.cutadapt:q} "
        "--raw-quality {input.raw:q} --output {output.table:q} --comparison-output {output.comparison:q}"

rule finalize_qc:
    input:
        tables=[qc + "/" + n for n in ("download_validation.tsv", "fastq_validation.tsv", "raw_read_quality.tsv", "trimmed_read_quality.tsv", "read_quality_comparison.tsv", "cutadapt_summary.tsv", "primer_detection_raw.tsv", "primer_detection_trimmed.tsv")],
        manifest=manifest_input, config="config/full_study.yaml", code=new_code,
    output:
        [qc + "/" + n for n in ("summary.json", "SUCCESS", "provenance.json", "input_checksums.json", "config_snapshot.yaml", "validation.json")],
    shell: "python scripts/study/qc_stage.py finalize --config {input.config:q}"

rule run_dada2:
    input:
        validated=qc + "/validation.json", raw=qc + "/raw_read_quality.tsv", trimmed=qc + "/trimmed_read_quality.tsv",
        manifest=manifest_input, config=ancient(scope["project_config"]), reads=lambda wc: reads(wc, "trimmed"),
        code=["scripts/dada2/run_pilot.R", "scripts/dada2/helpers.R"],
    output:
        [dada + "/" + n for n in ("SUCCESS", "read_tracking.tsv", "asv_counts.tsv", "asv_sequences.tsv", "asv_sequences.fasta", "summary.tsv", "sequence_length_distribution.tsv", "merge_diagnostics.tsv", "error_learning.tsv", "run_provenance.yaml", "input_checksums.tsv", "output_checksums.tsv", "config_snapshot.yaml", "sample_metadata.tsv", "software_versions.tsv", "session_info.txt", "warnings.txt")],
        figures=[dada + "/error_model_" + d + "." + ext for d in ("F", "R") for ext in ("pdf", "svg", "png")],
    params: output=dada, intermediate=scope["dada2_intermediate_dir"], reads=interim, settings=processing["dada2"], config_sha256=sha256_file(scope["project_config"]), manifest_sha256=manifest_hash
    threads: processing["dada2"]["threads"]
    resources: mem_mb=scope["resources"]["analysis_memory_mb"]
    log: scope["dada2_intermediate_dir"] + "/execution.log"
    shell:
        "mkdir -p {params.intermediate:q} && Rscript scripts/dada2/run_pilot.R --config {input.config:q} --manifest {input.manifest:q} "
        "--input-dir {params.reads:q} --output-dir {params.output:q} --intermediate-dir {params.intermediate:q} "
        "--raw-quality {input.raw:q} --trimmed-quality {input.trimmed:q} --threads {threads} > {log:q} 2>&1"

rule validate_dada2:
    input: products=rules.run_dada2.output, manifest=manifest_input, code=["scripts/dada2/validate_outputs.py", *analysis_code]
    output: dada + "/validation.json"
    params: directory=dada
    shell: "python scripts/dada2/validate_outputs.py --run-dir {params.directory:q} --check-inputs-root . --output {output:q}"

rule download_reference:
    input: config=ancient(scope["taxonomy_config"])
    output: protected(reference["path"])
    shell: "python scripts/taxonomy/download_reference.py --config {input.config:q} --database-only"

rule validate_reference:
    input: database=reference["path"], config=scope["taxonomy_config"], code="scripts/taxonomy/download_reference.py"
    output: reference["provenance_path"]
    params: config_sha256=sha256_file(scope["taxonomy_config"])
    shell: "python scripts/taxonomy/download_reference.py --config {input.config:q} --validate-only"

rule run_taxonomy:
    input:
        validated=dada + "/validation.json", config=scope["taxonomy_config"], reference=reference["path"], reference_provenance=reference["provenance_path"], manifest=manifest_input,
        data=[dada + "/" + n for n in ("SUCCESS", "asv_sequences.tsv", "asv_counts.tsv", "sample_metadata.tsv")],
        code=["scripts/taxonomy/run_taxonomy.py", "scripts/taxonomy/assign_taxonomy.R", "scripts/taxonomy/helpers.R", "scripts/dada2/helpers.R", "scripts/dada2/validate_outputs.py", "environment/conda-linux-64.lock", *analysis_code],
    output:
        [tax + "/" + n for n in ("SUCCESS", "provenance.json", "input_checksums.json", "taxonomy.tsv", "taxonomy_unfiltered.tsv", "taxonomy_bootstraps.tsv", "taxonomy_sensitivity.tsv", "taxonomy_screening.tsv", "assignment_coverage.tsv", "sample_coverage.tsv", "software_versions.tsv", "session_info.txt", "warnings.txt", "config_snapshot.yaml", "assignment_coverage.pdf", "assignment_coverage.svg", "assignment_coverage.png")],
    params: input_dir=dada, output_dir=tax, config_sha256=sha256_file(scope["taxonomy_config"])
    threads: 2
    resources: mem_mb=scope["resources"]["analysis_memory_mb"]
    shell: "python scripts/taxonomy/run_taxonomy.py --config {input.config:q} --input-dir {params.input_dir:q} --output-dir {params.output_dir:q} --threads {threads}"

rule validate_taxonomy:
    input: products=rules.run_taxonomy.output, manifest=manifest_input, code=["scripts/taxonomy/validate_taxonomy.py", "python/cacao_inventory/taxonomy_validation.py"]
    output: tax + "/validation.json"
    params: input_dir=dada, output_dir=tax
    shell: "python scripts/taxonomy/validate_taxonomy.py --input-dir {params.input_dir:q} --run-dir {params.output_dir:q} --output {output:q}"

rule bacterial_filter:
    input:
        validated=tax + "/validation.json", config=scope["diversity_config"], manifest=manifest_input,
        data=[dada + "/" + n for n in ("asv_sequences.tsv", "asv_counts.tsv", "sample_metadata.tsv")] + [tax + "/" + n for n in ("taxonomy.tsv", "taxonomy_bootstraps.tsv", "taxonomy_screening.tsv", "provenance.json")],
        code=["scripts/filtering/prepare_bacterial_table.py", "scripts/taxonomy/validate_taxonomy.py", "scripts/dada2/validate_outputs.py", "environment/conda-linux-64.lock", *analysis_code],
    output:
        [diversity["bacterial_dir"] + "/" + n for n in ("asv_counts.tsv", "asv_sequences.tsv", "taxonomy.tsv", "sample_metadata.tsv")],
        [diversity["filtering_dir"] + "/" + n for n in ("SUCCESS", "provenance.json", "input_checksums.json", "config_snapshot.yaml", "summary.json", "asv_filter_log.tsv", "sample_retention.tsv")],
    params: config_sha256=sha256_file(scope["diversity_config"])
    shell: "python scripts/filtering/prepare_bacterial_table.py --config {input.config:q}"

rule validate_filter:
    input: products=rules.bacterial_filter.output, config=scope["diversity_config"], manifest=manifest_input, code="scripts/filtering/validate_bacterial_table.py"
    output: diversity["filtering_dir"] + "/validation.json"
    shell: "python scripts/filtering/validate_bacterial_table.py --config {input.config:q} --output {output:q}"

rule describe_diversity:
    input:
        validated=diversity["filtering_dir"] + "/validation.json", config=scope["diversity_config"], manifest=manifest_input,
        data=[diversity["bacterial_dir"] + "/" + n for n in ("asv_counts.tsv", "asv_sequences.tsv", "taxonomy.tsv", "sample_metadata.tsv")],
        code=["scripts/diversity/run_diversity.py", "scripts/diversity/describe_pilot.R", "scripts/dada2/helpers.R", "scripts/filtering/validate_bacterial_table.py", "environment/conda-linux-64.lock", *analysis_code],
    output:
        [diversity["diversity_dir"] + "/" + n for n in ("SUCCESS", "provenance.json", "input_checksums.json", "config_snapshot.yaml", "alpha_diversity.tsv", "clr_coordinates.tsv", "aitchison_distances.tsv", "bray_curtis_distances.tsv", "pca_scores.tsv", "pca_variance.tsv", "software_versions.tsv", "session_info.txt", "warnings.txt")],
        [diversity["diversity_dir"] + "/" + n + "." + ext for n in ("alpha_diversity", "aitchison_pca") for ext in ("pdf", "svg", "png")],
    params: config_sha256=sha256_file(scope["diversity_config"])
    shell: "python scripts/diversity/run_diversity.py --config {input.config:q}"

rule validate_diversity:
    input: products=rules.describe_diversity.output, config=scope["diversity_config"], manifest=manifest_input, code=["scripts/diversity/validate_diversity.py", "python/cacao_inventory/diversity_validation.py"]
    output: diversity["diversity_dir"] + "/validation.json"
    shell: "python scripts/diversity/validate_diversity.py --config {input.config:q} --output {output:q}"
