#%%
# import argparse, os, glob
import pandas as pd

# import matplotlib.pyplot as plt
from datetime import datetime as dt
import time
import click
import numpy as np
import yaml
from itertools import chain


# %%
@click.command()
@click.option("--metadata-path", required=True)
@click.option("--mutation-summary-path", required=True)
@click.option("--output-path", required=True)
@click.option("--reference-accession", required=True)
def enrich(metadata_path, mutation_summary_path, output_path, reference_accession):  
    mutation_summary = pd.read_csv(
        mutation_summary_path,
        header=0,
        sep="\t",
    )

    mutation_summary = mutation_summary.assign(s1_dist=mutation_summary.nucleotide.apply(lambda x: len(str(x).split(","))))
    mutation_summary = mutation_summary.rename(columns={mutation_summary.columns[0]: "strainName"})

    meta = pd.read_csv(metadata_path, sep="\t")
    meta = meta.merge(mutation_summary[["strainName", "s1_dist"]], on="ncbiAcc", how="inner")

    def parse_isostring(string):
        try:
            split_string = list(map(int, string.split("-")))
        except:
            return (1970, 1, 1)
        length = len(split_string)
        if length == 1:
            date = (split_string[0], 7, 1)
        if length == 2:
            date = (split_string[0], split_string[1], 15)
        if length == 3:
            date = (split_string[0], split_string[1], split_string[2])
        return date

    def parse_ambiguous(string):
        try:
            split_string = list(map(str, string.split("-")))
        except:
            return ("XXXX", "XX", "XX")
        length = len(split_string)
        if length == 1:
            date = (split_string[0], "XX", "XX")
        if length == 2:
            date = (split_string[0], split_string[1], "XX")
        if length == 3:
            date = (split_string[0], split_string[1], split_string[2])
        return date

    # From https://stackoverflow.com/a/6451892/7483211
    def toYearFraction(date):
        def sinceEpoch(date):  # returns seconds since epoch
            return time.mktime(date.timetuple())

        s = sinceEpoch

        year = date.year
        startOfThisYear = dt(year=year, month=1, day=1)
        startOfNextYear = dt(year=year + 1, month=1, day=1)

        yearElapsed = s(date) - s(startOfThisYear)
        yearDuration = s(startOfNextYear) - s(startOfThisYear)
        fraction = yearElapsed / yearDuration

        return date.year + fraction

    meta["num_date"] = meta.date.apply(
        lambda x: toYearFraction(dt(*parse_isostring(x)))
    )
    meta["amb_date"] = meta.date.apply(lambda x: "-".join(parse_ambiguous(x)))

    # deduplicate
    keep = meta[meta.ncbiAcc == reference_accession]
    meta = pd.concat(
        [
            keep,
            meta[~meta.strainName.isin(keep.strainName.values)].drop_duplicates(
                ["strainName"]
            ),
        ]
    )

    meta = meta.drop(columns=['date']).rename(columns={'amb_date':'date'})

    meta.to_csv(output_path, sep="\t")


if __name__ == "__main__":
    enrich()
