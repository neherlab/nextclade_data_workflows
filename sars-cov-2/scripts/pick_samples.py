#%%
import pandas as pd
import numpy as np
import click

#%%
# Should input pangos with bad ones excluded
# Output pangos where there were insufficient samples
# Load problematic excludes
# des = pd.read_csv('pango_filtered.tsv', sep='\t', infer_datetime_format=True)
@click.command()
@click.option("--designations", default="open_pango_metadata.tsv")
@click.option("--counts", default="nr.tsv")
@click.option("--exclude", default="problematic_exclude.txt")
@click.option("--outfile", default="chosen_pango_strains.txt")
def pick_samples(designations, counts, exclude, outfile):
    des = pd.read_csv(
        designations,
        sep="\t",
        parse_dates=False,
        usecols=[
            "strain",
            "date",
            "region",
            "Nextstrain_clade",
            "pango_lineage",
            "clock_deviation",
        ],
    )
    # des
    #%%
    des = des.replace([np.inf, -np.inf], np.nan)
    des.dropna(subset=["clock_deviation"], inplace=True)
    #%%
    des.date = pd.to_datetime(des.date)
    # %%
    # Should still filter out exclusions from diagnostic.py
    # des.groupby('pango_lineage').clock_deviation.median().sort_values().to_csv('clock_deviation.tsv', sep='\t')
    # # %%
    # des.groupby('pango_lineage').date.quantile([0.05,0.5,0.95]).to_csv('clock_deviation_q.tsv', sep='\t')

    #%%
    excl = pd.read_csv(exclude, sep=" ", header=None, names=["strain"])
    #%%
    # Exclude strains in problematic_exclude.txt
    des = des[~des.strain.isin(excl.strain)]
    #%%
    lin = pd.read_csv(
        counts, sep="\t", names=["lineage", "counts"], skiprows=1, index_col=0
    )
    # lin
    #%%
    # lineage - counts, q5date
    lin = pd.merge(
        lin,
        des.groupby("pango_lineage").date.quantile(
            0.05, interpolation="nearest"
        ),
        left_index=True,
        right_index=True,
        how="outer",
    ).sort_values("date", ascending=True)

    #%%
    lin["count_with_recent"] = lin["counts"].copy()
    lin["count_with_recent"] = lin["count_with_recent"].fillna(0).astype(int)
    # lin
    #%%
    lin.loc[
        pd.Timestamp("today") - pd.Timedelta(days=360) < lin.date,
        "count_with_recent",
    ] = lin["count_with_recent"].apply(lambda x: max(1, x))
    # lin['count_with_recent'] = lin.count_with_recent.astype(int)
    # lin[-60:]
    #%%
    lineages = set()
    for i, g in des.groupby("pango_lineage"):
        picked = set()
        target = lin.loc[i, "count_with_recent"]
        if target > len(g.strain):
            print(
                f"{i} has {len(g.strain)} samples, but should sample {target}"
            )
            target = len(g.strain)
        # if target > 0:
        #     picked.add(g[g.date == lin.loc[i, "date"]].strain.iloc[0])
        #     if len(picked) != 1:
        #         print(f"{i} has {len(picked)} after quantile")
        # Commented out because synthetic used instead
        c = 0
        while len(picked) < target-1 and c < 10 * len(g.strain):
            picked.add(g.strain.sample(n=1).item())
            c += 1
        if len(picked) != target:
            print(f"{i} has {len(picked)}, but should be {target}")
        lineages |= picked
    # lineages

    # %%
    # Lineages lacking designated samples
    print(f"Missing in metadata.tsv:\n{lin[~lin.index.isin(des.pango_lineage.unique())].counts}")
    # %%
    pd.Series(list(lineages)).to_csv(
        outfile, sep="\t", index=False, header=False
    )


# %%
if __name__ == "__main__":
    pick_samples()
