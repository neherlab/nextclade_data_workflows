import pandas as pd
import json
import argparse
from Bio import SeqIO

def main(args):
    df = pd.read_csv(
        "https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_binding_vs_neut_escape/main/2022_Cao_convergent/convergent_RBD_evolution/bind_expr/bind_expr_BA2.csv"
    )
    sites = df["site"].unique()
    mutations = df["mutation"].unique()
    df.set_index(["site", "mutation"], inplace=True)
    wuhan = str(SeqIO.read(args.wuhan_spike, "fasta").seq)
    ref = str(SeqIO.read(args.reference_spike, "fasta").seq)

    result = {"name": "binding", "weight": 1.0, "locations": {}}
    for site in sites:
        wt = ref[site - 1]
        if wt == "-":
            wt = wuhan[site - 1]
        offset = df.loc[(site, wt), "bind_avg"]
        result["locations"][str(site - 1)] = {}
        for mutation in mutations:
            value = df.loc[(site, mutation), "bind_avg"]
            if pd.isna(value):
                value = 0
            else:
                value = round(value - offset, 5)
            result["locations"][str(site - 1)][mutation] = value


    final = {
      "aaRange": {
        "begin": 330,
        "end": 531
      },
      "description": "This column displays the predicted impact of RBD mutations on ACE2 binding, relative to a BA.2 baseline.\nA score is in log10 space. Hence a score of +1 should be interpreted as 10x higher affinity.\nThe score assumes that mutations do not interact. The underlying data looked at the impact of each mutation independently.\nThe score is calculated using the same data and logic as Bloom lab's ACE2 binding calculator (see https://doi.org/10.1101/2022.09.20.508745).",
      "cds": "S",
      "ignore": {
        "clades": [
          "outgroup"
        ]
      },
      "name": "ace2_binding",
      "nameFriendly": "ACE2 binding",
      "data": [
          result
      ]
    }

    json.dump(final, open(args.output, "w"), indent=2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate escape calculator parameter JSON"
    )
    parser.add_argument(
        "--wuhan-spike", required=True, help="Path to wuhan spike translation"
    )
    parser.add_argument(
        "--reference-spike", required=True, help="Path to reference spike translation"
    )
    parser.add_argument(
        "--output", required=True, help="Path for output escape JSON file"
    )

    args = parser.parse_args()
    main(args)

# %%

# df = pd.read_csv(
#     "https://raw.githubusercontent.com/jbloomlab/SARS2-mut-fitness/main/results/aa_fitness/aa_fitness.csv"
# )
# df = df[df.gene == "S"]
# sites = df["aa_site"].unique()
# df
# # %%
# ref = deepcopy(ba2)
# result = {}
# for site in df.aa_site.unique():
#     result[str(site - 1)] = {
#         "default": -3.0,
#     }
#     for mutation in df[df.aa_site == site].aa.unique():
#         value = df[(df.aa_site == site) & (df.aa == mutation)].fitness.values[0]
#         result[str(site - 1)][mutation] = value

# result = {
#     "aaRange": {"begin": 0, "end": 1272},
#     "data": [{"name": "S", "weight": 1.0, "locations": result}],
#     "description": "Bloom-Neher fitness score",
#     "cds": "S",
#     "name": "fitness",
#     "nameFriendly": "fitness",
# }
# result
# # %%
# FROM_PATH = "/Users/corneliusromer/code/nextclade_data/data/nextstrain/sars-cov-2/BA.2/pathogen.json"
# TO_PATH = "/Users/corneliusromer/code/nextclade_data/data/nextstrain/sars-cov-2/BA.2/pathogen.json"
# pathogen = json.load(open(TO_PATH))
# # Check if it contains any phenotypeData
# if "phenotypeData" not in pathogen:
#     # If not, add it
#     pathogen["phenotypeData"] = []
# # Check if it already contains phenotypeData item with name=="ace2_binding",
# if any([x["name"] == "fitness" for x in pathogen["phenotypeData"]]):
#     # If so, replace it
#     # Determine index of item with name=="ace2_binding"
#     index = [x["name"] for x in pathogen["phenotypeData"]].index("fitness")
#     pathogen["phenotypeData"][index]["data"] = [result]
# else:
#     # If not, append it
#     pathogen["phenotypeData"].append(result)

# json.dump(pathogen, open(TO_PATH, "w"), indent=2)


# # %%
