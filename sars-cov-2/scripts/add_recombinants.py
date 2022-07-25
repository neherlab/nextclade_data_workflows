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
    
    # Add common ancestor for all recombinants to prevent massive polytomy
    clade_cls = type(tree.root)
    rec_parent = clade_cls(name="rec_parent", branch_length=0.1)
    
    # Add each recombinant to root node
    for recombinant in recombinants:
        # Need large branch length to make sure recombinants ignored for ancestral reconstruction
        rec_parent.root.clades.append(Phylo.BaseTree.Clade(name=recombinant, branch_length=10.0))
    tree.root.clades.append(rec_parent)

    # Write tree to output
    Phylo.write(tree, output, "newick", format_branch_length='%1.8f', branch_length_only=True)

if __name__ == "__main__":
    typer.run(main)
