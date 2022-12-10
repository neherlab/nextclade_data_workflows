import argparse
import json
from Bio import SeqIO, Phylo


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="add root mutations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True)
    parser.add_argument('--reference', type=str, required=True)
    parser.add_argument('--translations', type=str,  nargs='+', required=True, help="amino acid alignment")
    parser.add_argument('--genes', type=str, nargs='+', required=True, help="gene names")
    parser.add_argument('--nuc-mutations', type=str, required=True)
    parser.add_argument('--aa-mutations', type=str, required=True)
    parser.add_argument('--output-nuc-mutations', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    parser.add_argument('--output-aa-mutations', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    # make sure we have a list of genes and translations
    genes = args.genes if type(args.genes)==list else [args.genes]
    translations = args.translations if type(args.translations)==list else [args.translations]

    # parse tree, get root name
    T = Phylo.read(args.tree, 'newick')
    root_name = T.root.name

    node_data = {}
    references ={}
    root_mutations = {'nuc':[], 'aa':{}}
    aa_references = {}

    with open(args.nuc_mutations, 'r') as fh:
        nuc_muts = json.load(fh)
    with open(args.aa_mutations, 'r') as fh:
        aa_muts = json.load(fh)

    root_seqs = aa_muts["nodes"][root_name]["aa_sequences"]

    # for each gene, get the root sequence and compare to reference
    for gene, translation in zip(genes, translations):
        seqs = {}
        for s in SeqIO.parse(translation, 'fasta'):
            seqs[s.id] = str(s.seq)

        root_mutations['aa'][gene] = []
        aa_references[gene] = seqs[args.reference]
        for pos, (ref, root) in enumerate(zip(seqs[args.reference], root_seqs[gene])):
            if ref!=root:
                root_mutations['aa'][gene].append(f"{ref}{pos+1}{root}")

    # repeat for nucleotide mutations
    for pos, (ref, root) in enumerate(zip(nuc_muts['nodes'][args.reference]['sequence'], nuc_muts['nodes'][root_name]['sequence'])):
        if ref!=root:
            root_mutations['nuc'].append(f"{ref}{pos+1}{root}")

    # attach differences as mutations to node data entry for root
    nuc_muts['nodes'][root_name]['muts'] = root_mutations['nuc']
    aa_muts['nodes'][root_name]['aa_muts'] = root_mutations['aa']
    # update the reference sequence of the root (which really is the parent of the root)
    aa_muts['reference'] = aa_references

    # output files
    with open(args.output_nuc_mutations, 'w') as fh:
        json.dump(nuc_muts, fh, indent=2)

    with open(args.output_aa_mutations, 'w') as fh:
        json.dump(aa_muts, fh, indent=2)

