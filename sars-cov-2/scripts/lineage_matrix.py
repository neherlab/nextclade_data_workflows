import typer


def main(
    ref: str = "reference.fasta",
    meta: str = "meta.tsv",
    out: str = "pango_matrix.npz",
):
    import numpy as np
    import pandas as pd
    from tqdm import tqdm

    df = pd.read_csv(meta, sep="\t")

    reference = open(ref, "r").read()

    def append_tails(start, end, missing):
        if pd.isna(missing):
            missing = ""
        if start == 1:
            missing = "1," + missing
        elif start > 1:
            missing = "1-" + str(int(start)) + "," + missing

        if end == len(reference) - 1:
            missing = missing + "," + str(len(reference))
        elif end < len(reference) - 1:
            missing = (
                missing + "," + str(int(end + 1)) + "-" + str(len(reference))
            )

        missing = missing.strip(",").replace(",,", ",")

        return missing

    df["missing"] = df.apply(
        lambda row: append_tails(
            row.alignmentStart, row.alignmentEnd, row.missing
        ),
        axis=1,
    )

    def char_to_int(char):
        if char == "A":
            return 0
        elif char == "C":
            return 1
        elif char == "G":
            return 2
        return 3


    ref_vec = np.zeros(6 * len(reference), dtype=np.int32)
    for i, c in enumerate(reference):
        ref_vec[i * 6 + char_to_int(c)] = 1

    def process_ranges(vec, ranges, index):
        if pd.isna(ranges):
            return vec
        for i in ranges.split(","):
            if i == "":
                continue
            split = i.split("-")
            if len(split) == 1:
                start = end = int(split[0])
            else:
                start = int(split[0])
                end = int(split[1])
            for j in range(start - 1, end):
                vec[6 * j : 6 * j + 5] = 0
                vec[6 * j + index] = 1
        return vec

    def process_missing(vec, missing):
        return process_ranges(vec, missing, 5)

    def process_deletions(vec, missing):
        return process_ranges(vec, missing, 4)

    def process_substitutions(vec, subs):
        if pd.isna(subs):
            return vec
        for sub in subs.split(","):
            pos = int(sub[1:-1]) - 1
            index = char_to_int(sub[-1])
            vec[6 * pos : 6 * pos + 5] = 0
            vec[6 * pos + index] = 1
        return vec

    def process_row(row):
        v = ref_vec.copy()
        process_missing(v, row["missing"])
        process_deletions(v, row["deletions"])
        process_substitutions(v, row["substitutions"])
        return v

    acc = np.zeros(6 * len(reference), dtype=np.int32)
    for _, row in df[df.pango_designated == "AY.4.2"].iterrows():
        acc += process_row(row)

    def process_group(group):
        acc = np.zeros(6 * len(reference), dtype=np.int32)
        for _, row in group.iterrows():
            acc += process_row(row)
        return acc

    result = {}
    for i, g in tqdm(df.groupby("pango_designated")):
        result[i] = process_group(g)

    np.savez_compressed(out, **result)


if __name__ == "__main__":
    typer.run(main)
