"""
Joins designations on strain names using canonical names
To canonicalize: convert alphabetical to uppercase
and throw out everything not alphanumeric nor slash, dash
"""

import concurrent
import concurrent.futures
import re

import click
import pandas as pd

ALPHA_DIGITS_REGEX = re.compile("[^A-Z0-9/-]")


def process_dataframe(filepath, sep, col_to_canonicalize, col_to_keep=None):
    """Reads, canonicalizes, and drops duplicates for one input DataFrame."""
    df = pd.read_csv(filepath, sep=sep)
    df["canonical"] = (
        df[col_to_canonicalize]
        .str.replace(ALPHA_DIGITS_REGEX, "", regex=True)
        .str.upper()
    )
    df.drop_duplicates(subset=["canonical"], keep="last", inplace=True)
    df.set_index("canonical", inplace=True)
    return df[[col_to_keep]]


@click.command()
@click.option(
    "--metadata-strainnames", default="pre-processed/metadata_strainnames.tsv"
)
@click.option("--pango-in", default="pre-processed/designations.csv")
@click.option("--pango-designations", default="pre-processed/pango_designations.csv")
@click.option(
    "--pango-designated-strains", default="pre-processed/pango_designated_strains.txt"
)
def format(
    metadata_strainnames, pango_in, pango_designations, pango_designated_strains
):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit processing tasks for both files
        future_meta = executor.submit(
            process_dataframe, metadata_strainnames, "\t", "strain", "strain"
        )
        future_pango = executor.submit(
            process_dataframe, pango_in, ",", "taxon", "lineage"
        )
        meta_strains = future_meta.result()
        pango_des = future_pango.result()

    joined = meta_strains.join(pango_des, how="inner")

    joined[["strain", "lineage"]].to_csv(pango_designations, index=False)
    joined[["strain"]].to_csv(pango_designated_strains, header=False, index=False)


if __name__ == "__main__":
    format()
