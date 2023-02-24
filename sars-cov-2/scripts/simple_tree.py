import typer
from Bio import SeqIO

def main(
    alignment: str = typer.Option(..., help="Alignment in fasta format"),
    tree: str = typer.Option(..., help="Tree in newick format"),
):
    # Read in alignment
    alignment = SeqIO.to_dict(SeqIO.parse(alignment, "fasta"))
    # Return strain names
    strains = list(alignment.keys())

    if len(strains) == 0:
        raise ValueError("No strains found in alignment")
    if len(strains) == 1:
        nwk = f"({strains[0]}:1.0);"
    if len(strains) == 2:
        nwk = f"({strains[0]}:1.0,{strains[1]}:1.0);"
    
    with open(tree, "w") as f:
        f.write(nwk)

if __name__ == "__main__":
    typer.run(main)
