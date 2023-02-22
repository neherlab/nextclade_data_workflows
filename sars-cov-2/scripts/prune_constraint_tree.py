"""
Prunes constraint tree to only include sequences in the alignment
Usage:
python3 scripts/prune_constraint_tree.py \
    --constraint-tree {input.constraint_tree} \
    --strains {input.strains} \
    --output {output.constraint_tree}
"""

import typer
from Bio import Phylo

def main(
    constraint_tree: str = typer.Option(..., help="Constraint tree"),
    strains: str = typer.Option(..., help="List of strains"),
    output: str = typer.Option(..., help="Output file"),
):
    strains = set(open(strains).read().splitlines())
    tree = Phylo.read(constraint_tree, "newick")
    # Remove all tips that are not in the strains list
    for clade in tree.get_terminals():
        if clade.name not in strains:
            tree.prune(clade)
    Phylo.write(tree, output, "newick")

if __name__ == "__main__":
    typer.run(main)