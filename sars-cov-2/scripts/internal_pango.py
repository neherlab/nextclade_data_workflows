#%%
import itertools

# import sys
# sys.path.insert(1, '/Users/cr/code/treetime')
# import treetime
import augur.utils
import Bio.Align
import click
import pandas as pd
from pango_aliasor.aliasor import Aliasor
from Bio import Phylo
from treetime import TreeAnc


@click.command()
@click.option("--designations", required=True, type=str)
@click.option("--synthetic", required=False, type=str)
@click.option("--tree", required=True, type=str)
@click.option("--alias", required=True, type=str)
@click.option("--output", required=True, type=str)
@click.option("--field-name", default="inferred_lineage")
def main(designations, tree, alias, synthetic, output, field_name):
    """
    Takes designation csv, nwk tree, and alias json
    Produces node.json with field-name (default: inferred_lineage)
    """
    #%%
    # Initialize aliasor
    aliasor = Aliasor(alias)
    #%%
    # Read in meta
    meta = pd.read_csv(designations, index_col=0)

    #%%
    # Read in meta
    meta = pd.read_csv(designations, index_col=0)
    #%%
    # Read in synthetic lineages
    synthetic = pd.read_csv(synthetic, header=None, names=["lineage"])
    synthetic.index = synthetic.lineage

    meta = pd.concat([meta, synthetic], axis=0)

    #%%
    # Read tree
    tree = Phylo.read(tree, format="newick")
    #%%
    # Get tip names
    tips = list(map(lambda x: x.name, tree.get_terminals()))
    #%%
    # Filter meta down to strains in tree
    meta = meta[meta.index.isin(tips)]
    #%%
    # Unalias lineage names
    meta["unaliased"] = meta["lineage"].apply(aliasor.uncompress)
    #%%
    # Get lineages present in tree
    lineages_unaliased = set(meta.unaliased.unique())
    print(f"{len(lineages_unaliased)} lineages present in tips")
    #%%
    def get_lineage_hierarchy(lineage):
        """
        Takes lineage and returns list including all parental lineages
        >>> get_lineage_hierarchy("A.B.C.D.E")
        ['A', 'A.B', 'A.B.C', 'A.B.C.D', 'A.B.C.D.E']
        """
        lineage_split = lineage.split(".")
        hierarchy = []
        if len(lineage_split) == 1:
            return lineage_split
        for i in range(len(lineage_split)):
            hierarchy.append(".".join(lineage_split[0 : i + 1]))
        return hierarchy

    #%%
    # Set of lineages, including intermediary not present directly
    characters = sorted(
        set(
            itertools.chain.from_iterable(
                map(get_lineage_hierarchy, lineages_unaliased)
            )
        )
    )
    #%%
    # Dictionary that maps character to position in pseudosequence
    mapping = {}
    for pos, character in enumerate(characters):
        mapping[character] = pos

    inverse_mapping = {v: k for k, v in mapping.items()}
    #%%
    def lineage_to_vector(lineage):
        """
        Turn lineage into binary vector
        Convention: A = 0, G = 1
        """
        character_list = get_lineage_hierarchy(lineage)
        vector = ["A"] * len(characters)
        for character in character_list:
            vector[mapping[character]] = "G"
        return "".join(vector)

    #%%
    # Create pseudo sequence vector for all tips
    meta["trait_vector"] = meta["unaliased"].apply(lineage_to_vector)
    #%%
    # Transform binary vector into BioSeq alignment
    alignment = Bio.Align.MultipleSeqAlignment(
        [
            Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(vector), id=strain, name=strain)
            for strain, vector in meta.trait_vector.iteritems()
        ]
    )
    #%%
    # Requires patched treetime >=0.8.6 to work with large number of undesignated tips
    # https://github.com/neherlab/treetime/issues/177
    tt = TreeAnc(tree=tree, aln=alignment, ignore_missing_alns=True)
    tt.infer_ancestral_sequences(
        method="ml", marginal=False, reconstruct_tip_states=True
    )
    #%%
    # Get reconstructed alignment
    aln = tt.get_tree_dict()
    #%%
    # Convert reconstructed alignment to lineage
    def seq_to_lineage_list(seq):
        """
        Takes a BioSeq object and returns a list of lineages
        >>> seq_to_lineage_list("GGAAG")
        ['B', 'B.1', 'B.1.3']
        """
        lineage_list = []
        for index, character in enumerate(seq):
            if character == "G":
                lineage_list.append(inverse_mapping[index])
        return lineage_list

    def lineage_list_to_lineage(lineage_list):
        return lineage_list[-1]

    #%%
    for rec in aln:
        meta.loc[rec.id, "reconstructed_full"] = str(
            seq_to_lineage_list(rec.seq)
        )
        meta.loc[rec.id, "reconstructed"] = lineage_list_to_lineage(
            seq_to_lineage_list(rec.seq)
        )
    #%%
    meta["realiased"] = meta["reconstructed"].apply(aliasor.compress)
    #%%
    augur.utils.write_json(
        {
            "nodes": meta["realiased"]
            .rename(field_name)
            .to_frame()  # Necessary to prevent slash escape
            .to_dict(orient="index")
        },
        output,
    )


if __name__ == "__main__":
    main()
    # For debugging
    # main(['--tree', 'builds/nextclade/tree.nwk', '--synthetic', 'builds/nextclade/chosen_synthetic_strains.txt', '--designations', 'pre-processed/pango_designations_nextstrain_names.csv', '--output', 'builds/nextclade/internal_pango.json', '--field-name', 'Nextclade_pango'])
