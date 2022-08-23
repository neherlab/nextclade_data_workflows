from collections import defaultdict
import json
from typing import Dict, Set, Tuple
import typer


def main(
    ref: str = "references/MN908947/reference.fasta",
    matrix: str = "pre-processed/pango_matrix.npz",
    alias: str = "",
    overwrites: str = "",
    out: str = "synthetic.fasta",
):
    import numpy as np
    from pango_aliasor.aliasor import Aliasor
    from Bio import SeqIO

    revs: defaultdict[str, list] = defaultdict(list)

    npzfile = np.load(matrix)

    reference = str(SeqIO.read(ref, "fasta").seq)

    overwrite_dict: Dict[str, Set[Tuple[int, str]]] = {}

    def char_to_int(char):
        if char == "A":
            return 0
        elif char == "C":
            return 1
        elif char == "G":
            return 2
        elif char == "T":
            return 3
        elif char == "-":
            return 4

    def int_to_char(integer):
        if integer == 0:
            return "A"
        elif integer == 1:
            return "C"
        elif integer == 2:
            return "G"
        elif integer == 3:
            return "T"
        return "-"

    ref_vec = np.zeros(len(reference), dtype=np.int8)
    for i, c in enumerate(reference):
        ref_vec[i] = char_to_int(c)

    def mutations(lineage):
        """
        Given a lineage, return a dict of mutations
        """
        muts = {}
        seq = npzfile[lineage]
        seqr = seq.reshape(-1, 6)[:, :5]
        seqx = np.argmax(seqr, axis=1)
        mut_pos = np.nonzero(seqx - ref_vec)[0]
        for i in mut_pos:
            if seqr[i, seqx[i]] > 0:
                muts[(i + 1, int_to_char(seqx[i]))] = (
                    seqr[i, seqx[i]] / seqr[i].sum()
                )
        return muts

    if alias != "":
        aliasor = Aliasor(alias)
    else:
        aliasor = Aliasor()

    def parent(lineage):
        """
        Returns parent lineage
        """
        uncompressed = aliasor.uncompress(lineage)
        lineage_split = uncompressed.split(".")
        if len(lineage_split) == 1:
            return None
        else:
            try:
                compressed = aliasor.compress(".".join(lineage_split[:-1]))
            except:
                print(
                    f"uncompressed: {uncompressed}, lineage_split: {lineage_split}, lineage_split[:-1]: {lineage_split[:-1]}"
                )
                raise
            return compressed

    def reversion_check(lineage, muts):
        """
        Check if mutations are reverted in lineage
        Takes in set of (pos, char) tuples
        Checks which have evidence for reversion
        """
        try:
            seq = npzfile[lineage]
        except KeyError:
            return {}

        seqr = seq.reshape(-1, 6)[:, :5]

        mut_dict = {}
        for (pos, char) in muts:
            mut_dict[(pos, char)] = (
                seqr[pos - 1, char_to_int(char)],
                seqr[pos - 1].sum(),
            )

        return mut_dict

    muts: dict[str, set] = {}

    def defining_mutations(lineage: str):
        """
        Returns list of mutations compared to parent lineage
        Method:
        1. Lookup parent lineage mutations
        2. Check for new mutations
        3. Double check reversions
        """
        parent_lineage = parent(lineage)
        if parent_lineage is None:
            lineage_muts = set()
        else:
            if parent_lineage not in muts:
                defining_mutations(parent_lineage)
            lineage_muts = set(muts[parent_lineage])

        for mut, (present, total) in reversion_check(
            lineage, lineage_muts
        ).items():
            if total > 3 and present / total < 0.2:
                revs[lineage].append(f"{mut}: ({present}, {total})")
                lineage_muts.remove(mut)

        try:
            lineage_mutations = mutations(lineage)
        except KeyError:
            lineage_mutations = {}

        for (pos, char), freq in lineage_mutations.items():
            if freq > 0.9:
                lineage_muts.add((pos, char))
        
        for (pos, char) in overwrite_dict.get(lineage, set()):
            # Remove any mutations at overwrite position
            new_lineage_muts = set(
                [
                    mut
                    for mut in lineage_muts
                    if mut[0] != pos
                ]
            )
            if char != reference[pos - 1]:
                new_lineage_muts.add((pos, char))
            
            # Printing to help understand overwrites
            if new_lineage_muts - lineage_muts:
                print(
                    f"Overwrote {lineage} with {pos}{char}, added: {new_lineage_muts - lineage_muts}"
                )
            if lineage_muts - new_lineage_muts:
                print(
                    f"Overwrote {lineage} with {pos}{char}, removed: {lineage_muts - new_lineage_muts}"
                )
            lineage_muts = new_lineage_muts

        muts[lineage] = lineage_muts
        return lineage_muts

    def clean_synthetic(lineage):
        """
        Given a lineage, return a synthetic sequence
        """
        template = list(reference)
        for (pos, char) in defining_mutations(lineage):
            template[pos - 1] = char
        return "".join(template)

    if overwrites != "":
        import pandas as pd

        overwrite_file = pd.read_csv(overwrites, sep="\t")
        for i, row in overwrite_file.iterrows():
            lineage = row["lineage"]
            if lineage not in overwrite_dict:
                overwrite_dict[lineage] = set()
            overwrite_dict[lineage].add((int(row["pos"]), row["char"]))

    with open(out, "w") as f:
        muts = {}
        for lineage in npzfile.files:
            f.write(f">{lineage}\n")
            f.write(clean_synthetic(lineage))
            f.write("\n")
    
    with open(out + ".revs", "w") as f:
        json.dump(revs, f)


if __name__ == "__main__":
    typer.run(main)
