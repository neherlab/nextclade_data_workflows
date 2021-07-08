#%%
# import argparse, os, glob
import pandas as pd
# import matplotlib.pyplot as plt

#%%
vic = pd.read_csv(
    "pre-processed/vic.mutation_summary.tsv",
    header=1,
    names=["ncbiAcc", "nucleotide"],
    sep="\t",
)
yam = pd.read_csv(
    "pre-processed/yam.mutation_summary.tsv",
    header=1,
    names=["ncbiAcc", "nucleotide"],
    sep="\t",
)
# %%
vic = vic.assign(vic_dist=vic.nucleotide.apply(lambda x: len(str(x).split(","))))
vic = vic.assign(yam_dist=yam.nucleotide.apply(lambda x: len(str(x).split(","))))

# %%
meta = pd.read_csv("pre-processed/metadata.tsv", sep="\t")
meta = meta.merge(vic[["ncbiAcc", "vic_dist", "yam_dist"]], on="ncbiAcc")
meta = meta.assign(vic_yam=lambda x: x.vic_dist - x.yam_dist)
meta["subtype"] = meta.apply(
    lambda row: "yam" if row.yam_dist < row.vic_dist else "vic", axis=1
)
# %%

# %%
def parse_isostring(string):
    split_string = list(map(int,string.split("-")))
    length = len(split_string)
    if(length==1):
        date = (split_string[0],7,1)
    if(length==2):
        date = (split_string[0],split_string[1],15)
    if(length==3):
        date = (split_string[0],split_string[1],split_string[2])
    return date
# %%
meta["year"] = meta.date.apply(lambda x: parse_isostring(x)[0])
#%%
meta.to_csv("pre-processed/metadata_enriched.tsv", sep="\t")