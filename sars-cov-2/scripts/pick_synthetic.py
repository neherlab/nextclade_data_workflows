import pandas as pd
import numpy as np
import click


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

    des = des.replace([np.inf, -np.inf], np.nan)
    des.dropna(subset=["clock_deviation"], inplace=True)

    des.date = pd.to_datetime(des.date)

    lin = pd.read_csv(
        counts, sep="\t", names=["lineage", "counts"], skiprows=1, index_col=0
    )

    lin = pd.merge(
        lin,
        des.groupby("pango_lineage").date.quantile(
            0.05, interpolation="nearest"
        ),
        left_index=True,
        right_index=True,
        how="outer",
    ).sort_values("date", ascending=True)

    lin["count_with_recent"] = lin["counts"].copy()
    lin["count_with_recent"] = lin["count_with_recent"].fillna(0).astype(int)

    lin.loc[
        pd.Timestamp("today") - pd.Timedelta(days=360) < lin.date,
        "count_with_recent",
    ] = lin["count_with_recent"].apply(lambda x: max(1, x))

    lineages = set()
    for i, g in des.groupby("pango_lineage"):
        target = lin.loc[i, "count_with_recent"]
        if target > 0:
            lineages.add(i)

    pd.Series(list(lineages)).to_csv(
        outfile, sep="\t", index=False, header=False
    )


if __name__ == "__main__":
    pick_samples()
