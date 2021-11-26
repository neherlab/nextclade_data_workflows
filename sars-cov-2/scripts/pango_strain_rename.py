import re

import click
import pandas as pd


@click.command()
@click.option("--metadata-strainnames", default="pre-processed/metadata_strainnames.tsv")
@click.option("--pango-in", default="pre-processed/pango_raw.csv")
@click.option("--pango-designations", default="pre-processed/pango_designations.csv")
@click.option("--pango-designated-strains", default="pre-processed/pango_designated_strains.txt")
def format(metadata_strainnames, pango_in, pango_designations, pango_designated_strains):
    meta_strains = pd.read_csv(metadata_strainnames, sep="\t")
    pango_des = pd.read_csv(pango_in, sep=",")

    alpha_digits_regex = re.compile("[^A-Z0-9/\-]")

    pango_des["canonical"] = pango_des.taxon.apply(
        lambda x: alpha_digits_regex.sub("", x.upper())
    )
    meta_strains["canonical"] = meta_strains.strain.apply(
        lambda x: alpha_digits_regex.sub("", x.upper())
    )

    pango_des.drop_duplicates(subset=["canonical"], keep=False, inplace=True)
    meta_strains.drop_duplicates(subset=["canonical"], keep=False, inplace=True)

    pango_des.set_index("canonical", inplace=True)
    meta_strains.set_index("canonical", inplace=True)

    pango_des = pango_des.join(meta_strains, on="canonical", how="left")

    pango_des.dropna(inplace=True)

    pango_des.to_csv(pango_designations, columns=["strain", "lineage"], index=False)
    pango_des.to_csv(
        pango_designated_strains, columns=["strain"], header=False, index=False
    )


if __name__ == "__main__":
    format()
