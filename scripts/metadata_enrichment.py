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

# #%%
# vic = pd.read_csv(
#     "pre-processed/vic.mutation_summary.tsv",
#     header=0,
#     sep="\t",
# )
# yam = pd.read_csv(
#     "pre-processed/yam.mutation_summary.tsv",
#     header=0,
#     sep="\t",
# )
# #%%
# vic = vic.assign(vic_dist=vic.nucleotide.apply(lambda x: len(str(x).split(","))))
# vic = vic.assign(yam_dist=yam.nucleotide.apply(lambda x: len(str(x).split(","))))
# vic = vic.rename(columns={vic.columns[0]: 'ncbiAcc'})
# # %%
# meta = pd.read_csv("pre-processed/metadata.tsv", sep="\t")
# meta = meta.merge(vic[["ncbiAcc", "vic_dist", "yam_dist"]], on="ncbiAcc")
# meta = meta.assign(vic_yam=lambda x: x.vic_dist - x.yam_dist)
# meta["subtype"] = meta.apply(
#     lambda row: "yam" if row.yam_dist < row.vic_dist else "vic", axis=1
# )
# # %%

# # %%
# def parse_isostring(string):
#     split_string = list(map(int,string.split("-")))
#     length = len(split_string)
#     if(length==1):
#         date = (split_string[0],7,1)
#     if(length==2):
#         date = (split_string[0],split_string[1],15)
#     if(length==3):
#         date = (split_string[0],split_string[1],split_string[2])
#     return date

# def parse_ambiguous(string):
#     split_string = list(map(str,string.split("-")))
#     length = len(split_string)
#     if(length==1):
#         date = (split_string[0],"XX","XX")
#     if(length==2):
#         date = (split_string[0],split_string[1],"XX")
#     if(length==3):
#         date = (split_string[0],split_string[1],split_string[2])
#     return date
# # %%
# # meta["year"] = meta.date.apply(lambda x: parse_isostring(x)[0])
# #%%
# # From https://stackoverflow.com/a/6451892/7483211
# def toYearFraction(date):
#     def sinceEpoch(date): # returns seconds since epoch
#         return time.mktime(date.timetuple())
#     s = sinceEpoch

#     year = date.year
#     startOfThisYear = dt(year=year, month=1, day=1)
#     startOfNextYear = dt(year=year+1, month=1, day=1)

#     yearElapsed = s(date) - s(startOfThisYear)
#     yearDuration = s(startOfNextYear) - s(startOfThisYear)
#     fraction = yearElapsed/yearDuration

#     return date.year + fraction
# #%%
# meta["num_date"] = meta.date.apply(lambda x: toYearFraction(dt(*parse_isostring(x))))
# meta["amb_date"] = meta.date.apply(lambda x: "-".join(parse_ambiguous(x)))
# #%%
# meta.to_csv("pre-processed/metadata_enriched.tsv", sep="\t")
# %%

# %%
@click.command()
@click.option("--flu-type", type=click.Choice(["A", "B"]), required=True)
def enrich(flu_type):
    if flu_type == "B":
        strain_1 = "vic"
        strain_2 = "yam"
    else:
        strain_1 = "h1"
        strain_2 = "h3"

    s1 = pd.read_csv(
        f"pre-processed/{strain_1}_4.mutation_summary.tsv",
        header=0,
        sep="\t",
    )
    s2 = pd.read_csv(
        f"pre-processed/{strain_2}_4.mutation_summary.tsv",
        header=0,
        sep="\t",
    )

    s1 = s1.assign(s1_dist=s1.nucleotide.apply(lambda x: len(str(x).split(","))))
    s2 = s2.assign(s2_dist=s2.nucleotide.apply(lambda x: len(str(x).split(","))))
    s1 = s1.rename(columns={s1.columns[0]: "ncbiAcc"})
    s2 = s2.rename(columns={s2.columns[0]: "ncbiAcc"})

    meta = pd.read_csv(f"pre-processed/metadata_{flu_type}_4.tsv", sep="\t")
    meta = meta.merge(s1[["ncbiAcc", "s1_dist"]], on="ncbiAcc", how="left")
    meta = meta.merge(s2[["ncbiAcc", "s2_dist"]], on="ncbiAcc", how="left").fillna(1000)
    meta = meta.assign(s1_s2=lambda x: x.s1_dist - x.s2_dist)
    meta["subtype"] = meta.apply(
        lambda row: strain_2 if row.s2_dist < row.s1_dist else strain_1, axis=1
    )

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
    with open("profiles/b/vic/builds.yaml", 'r') as stream:
        profile= yaml.safe_load(stream)
    r = profile['refine']['root']
    # roots = ["CY115151", "CY114381", "CY121680", "CY115183"]
    roots = [item for key in r.keys() for item in r[key].values() ]
    # print(roots)
    keep = meta[meta.ncbiAcc.isin(roots)]
    meta = pd.concat(
        [
            keep,
            meta[~meta.strainName.isin(keep.strainName.values)].drop_duplicates(
                ["strainName"]
            ),
        ]
    )

    meta.to_csv(f"pre-processed/metadata_enriched_{flu_type}.tsv", sep="\t")


if __name__ == "__main__":
    enrich()
