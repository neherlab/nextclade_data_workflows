#%%
import pandas as pd
import click


# %%
@click.command()
@click.option("--strain-names-path", required=True)
@click.option("--enriched-metadata-path", required=True)
@click.option("--output-path", required=True)
def map_strain_name_to_accessions(strain_names_path, enriched_metadata_path, output_path):  
    strain_names = pd.read_csv(
        strain_names_path,
        names=["strainName"],
        # header=None,
        sep="\t",
    )

    metadata = pd.read_csv(
        enriched_metadata_path,
        header=0,
        sep="\t",
    )

    strain_names = strain_names.merge(metadata, on="strainName", how="inner")
    strain_names.ncbiAcc.to_csv(output_path, sep="\t", index=False)

if __name__ == "__main__":
    map_strain_name_to_accessions()
