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
        (pos, char) -> freq, non-N count
        """
        muts = {}
        seq = npzfile[lineage]
        seqr = seq.reshape(-1, 6)[:, :5]
        seqx = np.argmax(seqr, axis=1)
        mut_pos = np.nonzero(seqx - ref_vec)[0]
        for i in mut_pos:
            if seqr[i, seqx[i]] > 0:
                muts[(i + 1, int_to_char(seqx[i]))] = (
                    seqr[i, seqx[i]] / seqr[i].sum(),
                    seqr[i].sum(),
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
        for pos, char in muts:
            mut_dict[(pos, char)] = (
                seqr[pos - 1, char_to_int(char)],
                seqr[pos - 1].sum(),
            )

        return mut_dict

    muts: dict[str, set] = {}

    def defining_mutations(lineage):
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

        for mut, (present, total) in reversion_check(lineage, lineage_muts).items():
            # Position 44 has artefactual flipping, require stricter thresholds
            if mut[0] == 44:
                if total > 15 and present / total < 0.1:
                    lineage_muts.remove(mut)
            elif total > 3 and present / total < 0.2:
                lineage_muts.remove(mut)

        try:
            lineage_mutations = mutations(lineage)
        except KeyError:
            lineage_mutations = {}

        for (pos, char), (freq, count) in lineage_mutations.items():
            if freq > 0.9 or freq >= 0.7 and count < 11:
                lineage_muts.add((pos, char))

        for pos, char in overwrite_dict.get(lineage, set()):
            # Remove any mutations at overwrite position
            new_lineage_muts = set([mut for mut in lineage_muts if mut[0] != pos])
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
        for pos, char in defining_mutations(lineage):
            template[pos - 1] = char
        return "".join(template)

    def parse_positions(pos_str):
        pos_str = str(pos_str)
        if "-" in pos_str:
            parts = pos_str.split("-", 1)
            if parts[0].isdigit() and parts[1].isdigit():
                return range(int(parts[0]), int(parts[1]) + 1)
        return [int(pos_str)]

    if overwrites != "":
        import pandas as pd

        overwrite_file = pd.read_csv(overwrites, sep="\t", dtype={"pos": str})

        wildcard_rows = []
        exact_rows = []
        for i, row in overwrite_file.iterrows():
            lineage = row["lineage"]
            if isinstance(lineage, str) and lineage.endswith("*"):
                wildcard_rows.append(row)
            elif isinstance(lineage, str) and lineage.strip():
                exact_rows.append(row)

        # Sort wildcards by depth (number of dots in uncompressed form)
        # General lineages first, more specific last — so specific overwrites general
        wildcard_rows.sort(
            key=lambda row: aliasor.uncompress(row["lineage"].rstrip("*")).count(".")
        )

        # Apply wildcards in topological order
        for row in wildcard_rows:
            base_lineage = row["lineage"].rstrip("*")
            base_uncompressed = aliasor.uncompress(base_lineage)
            positions = parse_positions(row["pos"])
            char = row["char"]
            for lin in npzfile.files:
                try:
                    lin_uncompressed = aliasor.uncompress(lin)
                except Exception:
                    continue
                if lin_uncompressed == base_uncompressed or lin_uncompressed.startswith(base_uncompressed + "."):
                    if lin not in overwrite_dict:
                        overwrite_dict[lin] = set()
                    for pos in positions:
                        overwrite_dict[lin] = set(
                            mut for mut in overwrite_dict[lin] if mut[0] != pos
                        )
                        overwrite_dict[lin].add((pos, char))

        # Apply exact entries last (highest priority)
        for row in exact_rows:
            lineage = row["lineage"]
            positions = parse_positions(row["pos"])
            char = row["char"]
            if lineage not in overwrite_dict:
                overwrite_dict[lineage] = set()
            for pos in positions:
                overwrite_dict[lineage] = set(
                    mut for mut in overwrite_dict[lineage] if mut[0] != pos
                )
                overwrite_dict[lineage].add((pos, char))

    with open(out, "w") as f:
        muts = {}
        for lineage in npzfile.files:
            f.write(f">{lineage}\n")
            f.write(clean_synthetic(lineage))
            f.write("\n")


if __name__ == "__main__":
    typer.run(main)
