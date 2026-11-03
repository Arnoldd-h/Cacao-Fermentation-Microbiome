"""Independent second-study DADA2; never rebuild baseline or terminal QC."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path("scripts/study_dada2").resolve()))
from execute import DADA_PRODUCTS, AUDIT_PRODUCTS, inputs, load_policy

policy_path = "config/prjna627078_dada2.yaml"
policy, terminal, scope, manifest, processing = load_policy(policy_path)

rule all:
    input: policy["audit_dir"] + "/validation.json"

rule independent_dada2:
    input: [str(path) for path in inputs(Path(policy_path), policy, terminal, scope, manifest)]
    output:
        [policy["output_dir"] + "/" + name for name in DADA_PRODUCTS] +
        [policy["audit_dir"] + "/" + name for name in
         [*AUDIT_PRODUCTS, "SUCCESS", "validation.json", "config_snapshot.yaml", "input_checksums.json", "provenance.json"]]
    threads: policy["threads"]
    resources: mem_mb=policy["memory_mb"]
    params: config=policy_path
    shell: "python scripts/study_dada2/execute.py execute --config {params.config:q}"
