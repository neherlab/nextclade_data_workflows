import typer


def main(
    metadata: str = "builds/nextclade/extracted_metadata.tsv",
    synthetic: str = "builds/nextclade/chosen_synthetic_strains.txt",
    outfile: str = "builds/nextclade/metadata.tsv",
):
    import pandas as pd

    metadata = pd.read_csv(metadata, sep="\t")

    strains = pd.read_csv(synthetic, sep="\t", header=None, names=["strain"])

    df = pd.concat([metadata, strains], axis=0)
    df.fillna("?", inplace=True)
    df.to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
