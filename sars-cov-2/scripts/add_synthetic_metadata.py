import typer


def main(
    synthetic: str = "builds/nextclade/chosen_synthetic_strains.txt",
    outfile: str = "builds/nextclade/metadata.tsv",
):
    import pandas as pd


    df = pd.read_csv(synthetic, sep="\t", header=None, names=["strain"])

    df.fillna("?", inplace=True)
    df.to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
