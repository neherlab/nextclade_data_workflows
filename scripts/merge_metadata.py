#%%
# import argparse, os, glob
import pandas as pd
# import matplotlib.pyplot as plt
from datetime import datetime as dt
import time
import click
import numpy as np
import yaml
#%%
@click.command()
@click.option("--flu-type", type=click.Choice(["A", "B"]), required=True)
@click.option("--segment", required=True)
def merge(flu_type, segment):
    a = pd.read_csv(
        f"pre-processed/metadata_{flu_type}_{segment}.tsv",
        header=0,
        sep="\t",
    )
    b = pd.read_csv(
        f"pre-processed/metadata_enriched_{flu_type}.tsv",
        header=0,
        sep="\t",
    )
    #%%
    meta = a.merge(b.iloc[:,-8:],on="strainName",how="inner")
    #%%
    with open("profiles/b/vic/builds.yaml", 'r') as stream:
            profile= yaml.safe_load(stream)
    r = profile['refine']['root']
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

    meta = meta.drop(columns=['strain_y','date']).rename(columns={'strain_x': 'strain','amb_date':'date'})

    meta.to_csv(f"pre-processed/metadata_enriched_{flu_type}_{segment}.tsv", sep="\t")
#%%

if __name__ == '__main__':
    merge()