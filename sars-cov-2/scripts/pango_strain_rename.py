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


def process_dataframe(filepath, sep, col_to_canonicalize):
    """Reads, canonicalizes, and drops duplicates for one input DataFrame."""
    df = pd.read_csv(filepath, sep=sep)
    df["canonical"] = (
        df[col_to_canonicalize]
        .astype(str)
        .str.upper()
        .str.replace(ALPHA_DIGITS_REGEX, "", regex=True)
    )
    df.drop_duplicates(subset=["canonical"], keep="first", inplace=True)

    # Determine columns to keep for merging
    cols_to_keep = ["canonical", col_to_canonicalize]
    if "lineage" in df.columns:  # Keep lineage only if it exists (for pango_des)
        cols_to_keep.append("lineage")
    return df[cols_to_keep]


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
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        # Submit processing tasks for both files
        future_meta = executor.submit(
            process_dataframe, metadata_strainnames, "\t", "strain"
        )
        future_pango = executor.submit(process_dataframe, pango_in, ",", "taxon")
        # Get results (processed DataFrames)
        meta_strains = future_meta.result()
        pango_des = future_pango.result()

    merged_data = pd.merge(
        pango_des[["lineage", "canonical"]],
        meta_strains[["strain", "canonical"]],
        on="canonical",
        how="inner",
    )

    merged_data[["strain", "lineage"]].to_csv(pango_designations, index=False)
    merged_data[["strain"]].to_csv(pango_designated_strains, header=False, index=False)


if __name__ == "__main__":
    format()
