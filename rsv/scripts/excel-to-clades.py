from typing import Annotated

import polars as pl
import typer


def main(
    excel: Annotated[str, typer.Option(help="Excel input")] = "/Users/corneliusromer/code/nextclade_data_workflows/rsv/profiles/consortium/amino-acid-genotypes.xlsx",
    sheet: Annotated[str, typer.Option(help="Sheet name")] = "AA_RSVA",
    outfile: Annotated[str, typer.Option(help="Clades.tsv output")] = "/Users/corneliusromer/code/nextclade_data_workflows/rsv/data/a/EPI_ISL_412866/clades_consortium_raw.tsv",
):
    df = pl.read_excel(excel, sheet_name=sheet)
    df = df.with_columns([
        pl.col("Signature Mutations").str.split(",").alias("mutations"),
    ])
    
    # Output into clades.tsv
    with open(outfile, "w") as f:
        f.write("clade\tgene\tsite\talt\n")
        for row in df.rows(named=True):
            if row["mutations"] is None:
                continue
            print(row)
            for mut in row["mutations"]:
                gene, rest = mut.split(":")
                site = rest[1:-1]
                alt = rest[-1]
                if gene == "G":
                    continue
                f.write(f"{row['RSV genotype']}\t{gene}\t{site}\t{alt}\n")

if __name__ == "__main__":
    typer.run(main)