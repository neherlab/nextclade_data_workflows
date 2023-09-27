import argparse
from pprint import pp

from Bio import Phylo


def get_terminals(clade):
    return frozenset(terminal.name for terminal in clade.get_terminals())


def get_splits_from_node_bio(clade):
    """
    Extract all splits from the given Bio.Phylo clade.
    """
    all_terminals = get_terminals(clade.root)
    splits = set()

    for internal_node in clade.find_clades(terminal=False):
        split_terminals = get_terminals(internal_node)
        splits.add(split_terminals)
        splits.add(all_terminals - split_terminals)

    return splits

def check_split_presence(constraint_split, splits):
    return constraint_split in splits
    # For each constraint split, check if one of the splits is a subset of one of the constraint splits


def check_violated_splits_bio(constraint_path, parsimony_path):
    constraint_tree = Phylo.read(constraint_path, "newick")
    parsimony_tree = Phylo.read(parsimony_path, "newick")

    # Prune all taxa that are not in the constraint tree
    taxa_to_prune = [
        terminal.name
        for terminal in parsimony_tree.get_terminals()
        if terminal.name not in (t.name for t in constraint_tree.get_terminals())
    ]


    for terminal in taxa_to_prune:
        parsimony_tree.prune(terminal)

    # all_taxa = frozenset(terminal.name for terminal in parsimony_tree.get_terminals())

    parsimony_splits = get_splits_from_node_bio(parsimony_tree.clade)
    constraint_splits = get_splits_from_node_bio(constraint_tree.clade)

    violated_splits = [
        split for split in constraint_splits if not check_split_presence(split, parsimony_splits)
    ]

    all_terminals = get_terminals(constraint_tree.clade)

    both_sided_violations = [
        (split, all_terminals - split) for split in violated_splits
    ]

    return both_sided_violations


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Check violated clades in unrooted trees using the split approach."
    )
    parser.add_argument(
        "constraint_tree_path",
        type=str,
        help="Path to the constraint tree in Newick format.",
    )
    parser.add_argument(
        "parsimony_tree_path",
        type=str,
        help="Path to the parsimony tree in Newick format.",
    )

    args = parser.parse_args()

    violations = check_violated_splits_bio(
        args.constraint_tree_path, args.parsimony_tree_path
    )
    pp(violations)
