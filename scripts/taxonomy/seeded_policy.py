"""Validate the optional registered adapter; absent means historical native RNG."""


def validate_native_rng(config: dict) -> None:
    policy = config["classification"].get("native_rng")
    if policy is not None and policy != {
            "adapter": "r_seeded_query_mt19937_v1", "dada2_version": "1.34.0",
            "upstream_git_commit": "5d453a416608f75f0bf50bcc74bbf0c1851e7f9f",
            "upstream_sha256": "1641cf17a322af65c31f6ad58fcf75e787ded90e9017960e5f9d3d2daa69a6a6"}:
        raise ValueError("Unregistered native taxonomy RNG adapter")
