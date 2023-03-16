"""
Add placement priors from nextclade.tsv to auspice.json tree

Usage:

python3 scripts/add_priors.py \
    --tree {input.tree} \
    --tsv {input.tsv} \
    --output {output.tree}
"""
import typer


def main(
    tree: str = "auspice/wuhan/auspice_raw.json",
    ndjson: str = "~/code/nextclade/out.tsv",
    output: str = "builds/wuhan/auspice_priors.json",
):
    """
    Add placement priors from nextclade.ndjson.zst to auspice.json tree
    """
    import json
    import polars as pl

    priors = (
        pl.scan_ndjson(ndjson, infer_schema_length=10000)
        .select(
            [
                pl.col("nearestNodes"),
                pl.col("nearestNodes")
                .arr.lengths()
                .pow(-1)
                .alias("1/nearestNodesListLength"),
            ]
        )
        .filter(pl.col("nearestNodes").arr.lengths() > 0)
        .explode("nearestNodes")
        .groupby("nearestNodes")
        .sum()
        .sort("1/nearestNodesListLength", descending=True)
        .with_columns(
            [
                (
                    pl.col("1/nearestNodesListLength")
                    / pl.col("1/nearestNodesListLength").sum()
                ).alias("priors"),
            ]
        )
        .with_columns(
            [
                pl.col("priors").log10().alias("log_priors"),
            ]
        )
        .collect()
    )
    print(priors)
    # Add priors to tree
    pd = {
        entry["nearestNodes"]: entry["log_priors"]
        for entry in priors.to_dicts()
    }
    node_data = pd

    auspice_json = json.load(open(tree, "r"))

    def attach_labels(n: dict[str, dict]):  # closure
        if "node_attrs" not in n:
            n["node_attrs"] = {}
        n["node_attrs"]["placement_prior"] = {"value": node_data.get(n["name"], -10)}

        if "children" in n:
            for c in n["children"]:
                attach_labels(c)

    attach_labels(auspice_json["tree"])

    auspice_json["meta"]["colorings"].append(
        {
            "key": "placement_prior",
            "title": "Placement Prior (log10)",
            "type": "continuous",
        }
    )

    json.dump(auspice_json, open(output, "w"), indent=2)


if __name__ == "__main__":
    typer.run(main)
