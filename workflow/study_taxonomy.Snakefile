"""Independent taxonomy and chronological repetition; consumes completed DADA2."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path("scripts/study_taxonomy").resolve()))
from execute import inputs, load_policy, products

policy_path = "config/prjna627078_taxonomy.yaml"
policy = load_policy(Path(policy_path))[0]

rule all:
    input: policy["audit_dir"] + "/validation.json"

rule independent_taxonomy:
    input: [str(path) for path in inputs(Path(policy_path))]
    output:
        [str(path) for path in products(policy)] +
        [policy["audit_dir"] + "/" + name for name in
         ["SUCCESS", "validation.json", "config_snapshot.yaml", "input_checksums.json", "provenance.json"]]
    threads: policy["threads"]
    resources: mem_mb=policy["memory_mb"]
    params: config=policy_path
    shell: "python scripts/study_taxonomy/execute.py execute --config {params.config:q}"
