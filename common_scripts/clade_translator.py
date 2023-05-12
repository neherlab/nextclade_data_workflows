#!/bin/env python3
import copy
import os

import numpy as np
from Bio import AlignIO, Seq, SeqIO, SeqRecord


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Annotate sequences using a genbank reference')
    parser.add_argument('--original-reference', help='Genbank sequence of the original reference')
    parser.add_argument('--new-reference', help='Genbank sequence of the new reference')
    parser.add_argument('--clades', help='clades tsv')
    parser.add_argument('--output-clades', required=True, type=str, help='Output file')
    return parser.parse_args()

def get_coordinate_map(ref, qry):
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as tmp_dir:
        SeqIO.write([SeqRecord.SeqRecord(id='ref', seq=Seq.Seq(ref)), SeqRecord.SeqRecord(id='qry', seq=Seq.Seq(qry))], f"{tmp_dir}/sequences.fasta", "fasta")

        os.system(f"mafft --quiet {tmp_dir}/sequences.fasta > {tmp_dir}/aligned.fasta")
        ref_aligned, qry_aligned = AlignIO.read(f"{tmp_dir}/aligned.fasta", "fasta")

    aln_to_ref = np.cumsum(np.array(ref_aligned.seq) !='-')
    aln_to_qry = np.cumsum(np.array(qry_aligned.seq) !='-')

    ref_to_qry = aln_to_qry[np.array(ref_aligned.seq) !='-'] - 1
    qry_to_ref = aln_to_ref[np.array(qry_aligned.seq) !='-'] - 1

    positions_mapping_to_end_of_query = np.where(ref_to_qry == len(qry)-1)[0]
    if len(positions_mapping_to_end_of_query) > 1:
        ref_to_qry[positions_mapping_to_end_of_query[1]:] = len(qry)

    positions_mapping_to_end_of_ref = np.where(qry_to_ref == len(ref)-1)[0]
    if len(positions_mapping_to_end_of_ref) > 1:
        qry_to_ref[positions_mapping_to_end_of_ref[1]:] = len(ref)

    return np.concatenate([ref_to_qry, [min(ref_to_qry[-1]+1, len(qry))]]), np.concatenate([qry_to_ref, [min(qry_to_ref[-1]+1, len(ref))]])


if __name__=="__main__":

    args = parse_args()

    orig_ref = SeqIO.read(args.original_reference, 'genbank')
    new_ref = SeqIO.read(args.new_reference, 'genbank')

    orig_features = {f.qualifiers['gene'][0]: f for f in orig_ref.features if 'gene' in f.qualifiers and f.type=='CDS'}
    new_features = {f.qualifiers['gene'][0]: f for f in new_ref.features if 'gene' in f.qualifiers and f.type=='CDS'}

    coord_maps = {'nuc': get_coordinate_map(orig_ref.seq, new_ref.seq)}
    for f in orig_features:
        try:
            coord_maps[f] = get_coordinate_map(orig_features[f].extract(orig_ref).translate().seq,
                                            new_features[f].extract(new_ref).translate().seq)
        except:
            print(f"Could not map {f}")


    with open(args.clades) as f:
        clades = [l.strip().split('\t') for l in f]

    with open(args.output_clades, 'w') as f:
        f.write("clade\tgene\tsite\talt\n")
        for clade in clades:
            if (len(clade) < 4) or clades[0][0]=='#':
                f.write('\t'.join(clade) + '\n')
                continue

            if clade[1] in coord_maps:
                try:
                    new_pos = max(0,coord_maps[clade[1]][0][int(clade[2])-1])+1
                    f.write('\t'.join([clade[0], clade[1], str(new_pos),clade[3]]) + '\n')
                except:
                    print(f"Could not map {clade}")
            
            if clade[1] == 'clade':
                f.write('\t'.join(clade) + '\n')

