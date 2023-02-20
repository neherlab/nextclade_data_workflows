import typer

import pandas as pd
from pango_aliasor.aliasor import Aliasor


def main(
    counts: str = "defaults/nr.tsv",
    lineages: str = "builds/nextclade/designated_lineages.txt",
    alias_file: str = "pre-processed/alias.json",
    outfile: str = "builds/nextclade/chosen_synthetic_strains.txt",
):
    aliasor = Aliasor(alias_file)

    def top_parent(lineage):
        while aliasor.parent(lineage) != "":
            lineage = aliasor.parent(lineage)
        return lineage

    lin = pd.read_csv(
        counts, sep="\t", names=["lineage", "counts"], skiprows=1, index_col=0
    )

    with open(lineages) as f:
        lineages = set(f.read().splitlines())

    keep = []
    for lineage in lineages:
        if (
            lineage in lin.index
            or aliasor.uncompress(lineage).startswith("B.1.1.529")
            or top_parent(lineage).startswith("X")
        ):
            keep.append(lineage)

    pd.Series(keep).to_csv(
        outfile, sep="\t", index=False, header=False
    )


if __name__ == "__main__":
    typer.run(main)
