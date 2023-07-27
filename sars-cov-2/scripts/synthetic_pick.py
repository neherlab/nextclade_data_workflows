import typer

import pandas as pd
from pango_aliasor.aliasor import Aliasor


def main(
    counts: str = "defaults/nr.tsv",
    lineages: str = "builds/nextclade/designated_lineages.txt",
    alias_file: str = "pre-processed/alias.json",
    build_name: str = "wuhan",  # other option is "21L"
    excluded_recombinants: str = "profiles/clades/21L/excluded_recombinants.txt",
    excludes: str = "profiles/clades/excludes.txt",
    outfile: str = "builds/nextclade/chosen_synthetic_strains.txt",
):
    aliasor = Aliasor(alias_file)

    def top_parent(lineage):
        while aliasor.parent(lineage) != "":
            lineage = aliasor.parent(lineage)
        return lineage

    excluded_recombinants = set(open(excluded_recombinants).read().splitlines())

    lin = pd.read_csv(
        counts, sep="\t", names=["lineage", "counts"], skiprows=1, index_col=0
    )

    with open(lineages) as f:
        lineages = set(f.read().splitlines())

    keep = []
    for lineage in lineages:
        if (
            (
                build_name == "wuhan" and lineage in lin.index
            )  # Include pre-Omicrons
            or aliasor.uncompress(lineage).startswith("B.1.1.529")
            or top_parent(lineage).startswith("X")
        ):
            if build_name == "21L":
                # Exclude BA.1 and BA.3s
                if aliasor.uncompress(lineage).startswith(
                    "B.1.1.529.1"
                ) or aliasor.uncompress(lineage).startswith("B.1.1.529.3"):
                    continue
                # Exclude non-BA.2 recombinants
                if top_parent(lineage) in excluded_recombinants:
                    continue
            if build_name == "22F":
                # Exclude all non-XBBs
                if not aliasor.uncompress(lineage).startswith("XBB"):
                    continue

            keep.append(lineage)
    
    # Add outgroup (BA.3)
    if build_name != "wuhan":
        keep.extend(["BA.3", "BA.1", "B.1.1", "B.1.617.2", "B.1", "B.1.1.7", "B", "A"])
    
    if build_name == "22F":
        keep.extend(["BA.2", "BA.4", "BA.5", "BQ.1", "BA.2.75", "CH.1.1", "BN.1", "BA.2.3.20", "BA.2.12.1"])
    
    # Remove excludes
    excludes = set(open(excludes).read().splitlines())
    keep = list(set(keep) - excludes)

    pd.Series(keep).to_csv(outfile, sep="\t", index=False, header=False)


if __name__ == "__main__":
    typer.run(main)
