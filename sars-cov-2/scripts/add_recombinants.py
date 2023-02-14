# This script is called like this:
# python add_recombinants.py \
# --tree {input.tree} \
# --recombinants {input.recombinants} \
# --root {input.root} \
# --output {output.tree} 2>&1 | tee {log}

import typer
import ipdb


def main(
    tree: str = "builds/nextclade/tree_raw.nwk",
    recombinants: str = "builds/nextclade/recombinants.txt",
    root: str = "Wuhan/Hu-1/2019",
    recombinant_trees: str = "",
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
    # Maybe also need to attach another reference here with 0 branch length
    # To prevent flip flopping

    def lookup_by_names(tree):
        names = {}
        for clade in tree.find_clades():
            if clade.name:
                if clade.name in names:
                    raise ValueError("Duplicate key: %s" % clade.name)
                names[clade.name] = clade
        return names
    
    # Create recombinant tree based on pango structure, start with top level, then add first level etc.
    # Attach each clade to parent, parent can be calculated by removing one `.`
    def parent(child: str) -> str:
        split = child.split(".")
        if len(split) == 1:
            return None
        return ".".join(split[:-1])
    
    def attach(tip: str):
        # Check if tip already in tree
        for clade in rec_parent.root.find_clades():
            if clade.name == tip:
                return
        
        if parent(tip) is None:
            rec_parent.root.clades.append(Phylo.BaseTree.Clade(name=tip, branch_length=10.0))
        else:
            # ipdb.set_trace()
            # check if internal parent node exists
            # if exists, attach there
            # if doesn't exist, create
            # before attaching, check if internal node with this name exists

            # check if parent internal node exists
            if f"internal_{parent(tip)}" not in lookup_by_names(rec_parent):
                # check if parent tip exists
                if parent(tip) not in lookup_by_names(rec_parent):
                    attach(parent(tip))
                    # rename parent tip to internal node, and attach parent tip to internal node
                lookup_by_names(rec_parent)[parent(tip)].name = f"internal_{parent(tip)}"
                lookup_by_names(rec_parent)[f"internal_{parent(tip)}"].clades.append(Phylo.BaseTree.Clade(name=parent(tip), branch_length=0))
            lookup_by_names(rec_parent)[f"internal_{parent(tip)}"].clades.append(Phylo.BaseTree.Clade(name=tip, branch_length=1/30000))
    
    # Attach recombinant trees
    for recombinant_tree in recombinant_trees.split(","):
        # import ipdb; ipdb.set_trace()
        # Get recombinant name (top parent)
        recombinant_tree = Phylo.read(recombinant_tree, "newick")
        rec_name = [r.name for r in recombinant_tree.get_terminals()][0]
        while parent(rec_name) is not None:
            rec_name = parent(rec_name)
        print(rec_name)
        recombinant_tree.root_with_outgroup(rec_name)
        recombinant_tree.root.name = f"internal_{rec_name}"
        recombinant_tree.root.branch_length=10.0
        rec_parent.root.clades.append(recombinant_tree.root)
        # lookup_by_names(rec_parent)[f"internal_{rec_name}"].clades.append(recombinant_tree.root)

    # Add each recombinant to root node
    for recombinant in recombinants:
        # Call attach child to parent recursively 
        attach(recombinant)

    # Phylo.draw_ascii(rec_parent)
        
    tree.root.clades.append(rec_parent)

    # Write tree to output
    Phylo.write(tree, output, "newick", format_branch_length='%1.8f', branch_length_only=True)

if __name__ == "__main__":
    typer.run(main)
