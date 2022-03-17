# This script is called like this:
# python add_recombinants.py \
# --tree {input.tree} \
# --recombinants {input.recombinants} \
# --root {input.root} \
# --output {output.tree} 2>&1 | tee {log}

import typer


def main(
    tree: str = "builds/nextclade/tree_raw.nwk",
    recombinants: str = "builds/nextclade/recombinants.txt",
    root: str = "Wuhan/Hu-1/2019",
    output: str = "builds/nextclade/tree_with_recombinants.nwk",
):
    from Bio import Phylo 

    # Load tree
    tree = Phylo.read(tree, "newick")

    # Root tree on reference
    tree.root_with_outgroup(root)

    # Load recombinants from txt into a list
    with open(recombinants, "r") as f:
        recombinants = f.read().splitlines()
    
    # Add each recombinant to root node
    for recombinant in recombinants:
        try:
            tree.prune(target=recombinant)
        except ValueError:
            pass
        tree.root.clades.append(Phylo.BaseTree.Clade(name=recombinant, branch_length=0.0001))

    # Write tree to output
    Phylo.write(tree, output, "newick", format_branch_length='%1.8f', branch_length_only=True)

if __name__ == "__main__":
    typer.run(main)
