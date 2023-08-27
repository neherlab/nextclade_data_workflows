from Bio import Phylo
import pandas as pd
import argparse
from collections import defaultdict
import json

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Assign clades to a tree')
    parser.add_argument('--tree', type=str, help='Newick tree file')
    parser.add_argument('--metadata', type=str, help='metadata file with subtype information')
    parser.add_argument('--clade-key', type=str, help='column name in metadata file with clade information')
    parser.add_argument('--id-field', type=str, help='column name in metadata file with node id')
    parser.add_argument('--root-clade', type=str, default='unassigned', help='clade for the root')
    parser.add_argument('--output', type=str, help='output file')
    args = parser.parse_args()

    tree = Phylo.read(args.tree, 'newick')
    metadata = pd.read_csv(args.metadata, sep='\t')

    node_to_subtype = {x[args.id_field]: x[args.clade_key] for _,x in metadata.iterrows()}
    nodes_by_subtype = defaultdict(list)

    for node in tree.get_terminals():
        nodes_by_subtype[node_to_subtype.get(node.name, None)].append(node)
    for subtype, nodes in nodes_by_subtype.items():
        print("Subtype {} has {} nodes".format(subtype, len(nodes)))

    clade_demarcations = {}
    for subtype, nodes in nodes_by_subtype.items():
        if subtype is None:
            continue
        clade = tree.is_monophyletic(nodes)
        if clade:
            clade_demarcations[subtype] = clade
        else:
            clade_demarcations[subtype] = tree.common_ancestor(nodes)
            print("Subtype {} is not monophyletic".format(subtype))
        clade_demarcations[subtype].clade_label = subtype

    for node in tree.find_clades(order = 'preorder'):
        node.clade = args.root_clade

    for node in tree.find_clades(order = 'preorder'):
        if hasattr(node, "clade_label"):
            node.clade = node.clade_label
            for child in node.find_clades():
                child.clade = node.clade_label

    node_data = {'nodes': {}}
    for node in tree.find_clades(order = 'preorder'):
        if node.is_terminal() and node.name in node_to_subtype:
            if node.clade != node_to_subtype[node.name]:
                print("Node {} has clade {} but subtype {}".format(node.name, node.clade, node_to_subtype[node.name]))
                node.clade = node_to_subtype[node.name]
        datum = {'clade_membership': node.clade}
        if hasattr(node, 'clade_label'):
            datum['clade_annotation'] = node.clade_label
        node_data['nodes'][node.name] = datum


    with open(args.output, 'w') as f:
        json.dump(node_data, f, indent=2)