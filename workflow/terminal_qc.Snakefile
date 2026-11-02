"""Second-pass terminal QC; consumes but never rebuilds the validated baseline."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path("scripts/terminal_qc").resolve()))
from run import HELPERS, REPORTS, load_policy, read_paths

configfile: "config/prjna627078_terminal_qc.yaml"
policy_path = "config/prjna627078_terminal_qc.yaml"
policy, scope, manifest, processing = load_policy(policy_path)
runs = {row["run_accession"]: row for row in manifest}
interim = policy["interim_dir"]
intermediate = policy["intermediate_dir"]
output_dir = policy["output_dir"]
code = [str(p) for p in Path("scripts/terminal_qc").glob("*.py")] + ["python/cacao_inventory/" + name + ".py" for name in HELPERS]

rule all:
    input: output_dir + "/validation.json"

rule terminal_trim:
    input:
        reads=lambda w: [str(p.relative_to(Path.cwd())) for p in read_paths(policy, scope, runs[w.run])[0]],
        code=code, config=policy_path, scope=policy["source_scope"], project=scope["project_config"],
        manifest=scope["manifest"], validation=scope["qc_dir"] + "/validation.json",
        provenance=scope["qc_dir"] + "/provenance.json", environment="environment/conda-linux-64.lock",
        workflow="workflow/terminal_qc.Snakefile"
    output:
        reads=expand(interim + "/{{run}}/{{run}}_{mate}.fastq.gz", mate=[1, 2]),
        audit=intermediate + "/{run}/audit.json",
        reports=expand(intermediate + "/{{run}}/{name}", name=REPORTS)
    resources: mem_mb=300
    shell: "python scripts/terminal_qc/run.py trim --config {input.config:q} --run-accession {wildcards.run:q}"

rule terminal_fastqc:
    input:
        reads=expand(interim + "/{{run}}/{{run}}_{mate}.fastq.gz", mate=[1, 2]),
        config=policy_path, environment="environment/conda-linux-64.lock"
    output:
        expand(intermediate + "/fastqc/{{run}}_{mate}_fastqc.{ext}", mate=[1, 2], ext=["zip", "html"])
    threads: policy["fastqc_threads"]
    resources: mem_mb=policy["fastqc_memory_mb"]
    params: directory=intermediate + "/fastqc"
    shell: "mkdir -p {params.directory:q} && fastqc --threads {threads} --outdir {params.directory:q} {input.reads:q}"

rule terminal_finalize:
    input:
        audits=expand(intermediate + "/{run}/{name}", run=list(runs), name=["audit.json", *REPORTS]),
        quality=expand(intermediate + "/fastqc/{run}_{mate}_fastqc.zip", run=list(runs), mate=[1, 2]),
        baseline=expand(scope["qc_dir"] + "/fastqc_trimmed/{run}_{mate}_fastqc.zip", run=list(runs), mate=[1, 2]),
        config=policy_path, code=code, workflow="workflow/terminal_qc.Snakefile"
    output:
        [output_dir + "/" + name for name in ("run_summary.tsv", "read_length_counts.tsv", "quality_coverage.tsv", "summary.json",
         "config_snapshot.yaml", "input_checksums.json", "provenance.json", "SUCCESS", "validation.json")]
    resources: mem_mb=600
    shell: "python scripts/terminal_qc/run.py finalize --config {input.config:q}"
