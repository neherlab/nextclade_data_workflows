import argparse

import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", type=str, required=True)
    parser.add_argument("--designations", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t")
    designations = pd.read_csv(args.designations, sep=",", names=["strain", "designation"])

    metadata["pango_designated"] = metadata.join(designations.set_index("strain"), on="strain")[
        "designation"
    ]
    metadata["pango_lineage"] = metadata["pango_designated"].fillna(metadata["pango_lineage"])
    metadata["recombinant"] = metadata["pango_lineage"].str.startswith("X")
    metadata.to_csv(args.output, sep="\t", index=False)
