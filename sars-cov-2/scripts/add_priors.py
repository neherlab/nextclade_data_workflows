"""
Add placement priors from nextclade.tsv to auspice.json tree

Usage:

python3 scripts/add_priors.py \
    --tree {input.tree} \
    --tsv {input.tsv} \
    --output {output.tree}
"""
import importlib.metadata as metadata
import json

import packaging.version as version
import polars as pl
import typer

app = typer.Typer(pretty_exceptions_enable=False)


@app.command()
def main(
    tree: str = "auspice/wuhan/auspice_raw.json",
    ndjson: str = "~/code/nextclade/out.tsv",
    output: str = "builds/wuhan/auspice_priors.json",
):
    """
    Add placement priors from nextclade.ndjson.zst to auspice.json tree
    """

    
    polars_version = version.parse(metadata.version("polars"))

    if polars_version < version.parse("0.18.0"):
        # rename `arr` to `list` for compatibility with polars < 0.18.0
        pl.Expr.list = pl.Expr.arr

    priors = (
        pl.scan_ndjson(ndjson, infer_schema_length=10000)
        .select(
            [
                pl.col("nearestNodes"),
                pl.col("nearestNodes")
                .list.len().cast(pl.Float32)
                .pow(-1).inspect()
                .alias("1/nearestNodesListLength"),
            ]
        )
        .inspect()
        .filter(pl.col("nearestNodes").list.len() > 0)
        .explode("nearestNodes")
        .group_by("nearestNodes")
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

    def attach_labels(n):  # closure
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
    app()