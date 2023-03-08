"""
Add placement priors from nextclade.tsv to auspice.json tree

Usage:

python3 scripts/add_priors.py \
    --tree {input.tree} \
    --tsv {input.tsv} \
    --output {output.tree}
"""
#%%
import typer

import polars as pl


def main(
    tree: str = "auspice/wuhan/auspice_raw.json",
    tsv: str = "~/code/nextclade/out.tsv",
    output: str = "builds/wuhan/auspice_priors.json",
):
    """
    Add placement priors from nextclade.tsv to auspice.json tree
    """
    # tree = pl.read_json(tree)

#%%
tsv: str = "~/code/nextclade/out.tsv"

pl.read_csv(tsv, sep="\t", infer_schema_length=10000).head(10)
#%%
priors = (
    pl.scan_csv(tsv, sep="\t", infer_schema_length=10000)
    .select(
        [
            pl.col("nearestNodes").str.split(";").alias("nearestNodes"),
            pl.col("nearestNodes")
            .str.split(";")
            .arr.lengths()
            .pow(-1)
            .alias("1/nearestNodesListLength"),
        ]
    )
    .explode("nearestNodes")
    .groupby("nearestNodes")
    .sum()
    .sort("1/nearestNodesListLength", descending=True)
    .with_columns(
        [
            pl.col("1/nearestNodesListLength").log10().alias("log_priors"),
        ]
    )
    .collect()
)
#%%
# Add priors to tree
pd = {entry["nearestNodes"]: entry["log_priors"] for entry in priors.to_dicts()}
node_data = pd

#

#%%
import json
auspice_json = json.load(open("auspice/wuhan/auspice_raw.json"))

#%%
def attach_labels(n: dict[str,dict]):  # closure
    if n["name"] in node_data:
        if "node_attrs" not in n:
            n["node_attrs"] = {}
        n["node_attrs"]["placement_prior"] = {"value": node_data[n["name"]]}

    if "children" in n:
        for c in n["children"]:
            attach_labels(c)

attach_labels(auspice_json["tree"])

with open("prior.json", "w") as f:
    json.dump(auspice_json, f, indent=2)

#%%
if __name__ == "__main__":
    typer.run(main)

#%%
import os
# print working directory
print(os.getcwd())
# %%
os.getcwd()

# %%
