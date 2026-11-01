"""Within-study inference after the separately completed full-study workflow.

Run workflow/study.Snakefile first if validated bacterial inputs are absent.
This workflow never rebuilds or modifies upstream scientific artifacts.
"""

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path("python").resolve()))
from cacao_inventory.config import load_json_yaml

configfile: "config/longitudinal.yaml"

source_config = load_json_yaml(config["source_diversity_config"])
source = source_config["paths"]
output_dir = config["output_dir"]
longitudinal_code = [str(path) for path in Path("scripts/longitudinal").glob("*") if path.suffix in (".py", ".R")]
shared_code = ["python/cacao_inventory/" + name + ".py" for name in ("analysis_artifacts", "bacterial_filter", "config", "io", "provenance", "taxonomy_reference", "taxonomy_validation")]
products = ["observations.tsv", "sample_membership.tsv", "batch_models.tsv", "slopes.tsv", "intercepts.tsv", "global_tests.tsv", "direction_cosines.tsv",
            "leave_one_batch_out.tsv", "sign_statistics.tsv", "trajectory_plot_data.tsv", "summary.json", "independent_validation.tsv", "software_versions.tsv",
            "session_info.txt", "community_trajectories.pdf", "community_trajectories.svg", "community_trajectories.png", "config_snapshot.yaml", "input_checksums.json", "provenance.json", "SUCCESS", "validation.json"]


rule all:
    input: output_dir + "/validation.json"


rule run_longitudinal:
    input:
        config="config/longitudinal.yaml", source_config=config["source_diversity_config"], project="config/config.yaml",
        source=[source["bacterial_dir"] + "/" + name for name in ("asv_counts.tsv", "sample_metadata.tsv", "taxonomy.tsv")],
        validated=[source["filtering_dir"] + "/" + name for name in ("SUCCESS", "validation.json", "provenance.json")],
        code=[*longitudinal_code, *shared_code, "workflow/longitudinal.Snakefile", "scripts/filtering/validate_bacterial_table.py", "scripts/taxonomy/validate_taxonomy.py", "scripts/taxonomy/seeded_policy.py"],
        environment="environment/conda-linux-64.lock",
    output: [output_dir + "/" + name for name in products]
    params: config_sha256=hashlib.sha256(Path("config/longitudinal.yaml").read_bytes()).hexdigest()
    resources: mem_mb=1000
    shell: "python scripts/longitudinal/run_longitudinal.py --config {input.config:q}"
