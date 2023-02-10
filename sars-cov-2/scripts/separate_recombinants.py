"""
Usage:
python3 scripts/separate-recombinants.py \
    --alignment {input.alignment} \
    --output-without-recombinants {output.without_recombinants} \
    --alias-json {input.alias_json} \
    --tree-recombinants {params.tree_recombinants}

Filters the alignment file into:
- File without any recombinants
- One file per recombinant passed via --tree-recombinants
"""

import typer
from Bio import SeqIO
from pango_aliasor.aliasor import Aliasor

def main(
    alignment: str = typer.Option(..., help="Alignment file"),
    output_without_recombinants: str = typer.Option(..., help="Output file without recombinants"),
    alias_json: str = typer.Option(..., help="Path to JSON file with aliases"),
    tree_recombinants: str = typer.Option(..., help="Recombinants for which to produce separate alignment files"),
    recombinants: str = typer.Option(..., help="List of recombinants"),
):
    aliasor = Aliasor(alias_json)

    def top_level(lineage):
        if "." in lineage:
            return aliasor.uncompress(lineage).split(".")[0]
        return lineage

    # Read alignment
    with open(alignment, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    
    # Produce output without recombinants
    with open(output_without_recombinants, "w") as handle:
        for record in records:
            if not top_level(record.id).startswith("X"):
                SeqIO.write(record, handle, "fasta")
    
    # Produce list of recombinants
    with open(recombinants, "w") as handle:
        for record in records:
            if top_level(record.id).startswith("X"):
                handle.write(f"{top_level(record.id)}\n")
    
    # Produce separate alignment files for tree-recombinants
    for recombinant in tree_recombinants.split(","):
        with open(f"builds/nextclade/masked_recombinant_{recombinant}.fasta", "w") as handle:
            for record in records:
                if top_level(record.id) == recombinant:
                    SeqIO.write(record, handle, "fasta")
    

if __name__ == "__main__":
    typer.run(main)