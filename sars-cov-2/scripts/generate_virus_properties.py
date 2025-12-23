#%%
import json
from collections import defaultdict

import pandas as pd
from pango_aliasor.aliasor import Aliasor
from tqdm import tqdm

#%%
# Defines minimum proportions and counts for relevant mutations to keep
MIN_PROPORTION = 0.2
# MIN_COUNT = 100000

clades_with_high_proportion_threshold = [
    "21A",
    "21M",
    "21H",
    "21G",
    "21F",
    "21E",
    "21C",
    "21B",
    # "20I",
    "20H",
    "20D",
    "20C",
    "20B",
    # "20A",
    "19B",
    "19A",
]

HIGH_THRESHOLD_PROPORTION = 0.5
#%%
# aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.zst .
df = pd.read_csv(
    "metadata.tsv.zst",
    # "s3://nextstrain-ncov-private/metadata.tsv.zst",
    sep="\t",
    usecols=["clade_nextstrain", "Nextclade_pango", "substitutions"],
    parse_dates=False,
    dtype={"clade_nextstrain": str, "Nextclade_pango": str, "substitutions": str},
).dropna()
df.rename(columns={"clade_nextstrain": "Nextstrain_clade"}, inplace=True)


#%%
# Make unaliased pango column

aliasor = Aliasor()
df["unaliased"] = df["Nextclade_pango"].apply(aliasor.uncompress)
# %%
def accumulate_mutations(acc: defaultdict(int), row) -> defaultdict(int):
    try:
        for mutation in str(row).split(","):
            acc[mutation] += 1
    except:
        print(row)
        raise
    return acc


def aggregate_mutations(series) -> defaultdict(int):
    mutations = defaultdict(int)
    for row in tqdm(series):
        mutations = accumulate_mutations(mutations, row)
    return mutations

#%%
# Overwrite with new clade name when new clade not yet in data
# Set to 23C if it starts with XBB.1.9
# df.loc[df["unaliased"].str.startswith("B.1.1.529.2.75.3.4.1.1.1.1"),"Nextstrain_clade"] = "23C"
# df.loc[df["unaliased"].str.startswith("XBB.1.9"),"Nextstrain_clade"] = "23D"
# df.loc[df["unaliased"].str.startswith("B.1.1.529.2.86.1.1.16.1.7"),"Nextstrain_clade"] = "24H"
# df.loc[df["unaliased"].str.startswith("B.1.1.529.2.86.1.1.49.1.1.1.1.1"),"Nextstrain_clade"] = "24I"
# New clades 25D-25I (25B and 25C are now in ingest data)
# Use unaliased versions
df.loc[df["unaliased"].str.startswith("B.1.1.529.2.86.1.1.11.1.3.1.1.10.2.1"),"Nextstrain_clade"] = "25D"  # MC.10.2.1
df.loc[df["unaliased"].str.startswith("B.1.1.529.2.86.1.1.16.1.7.9.1.1"),"Nextstrain_clade"] = "25E"  # PY.1
df.loc[df["unaliased"].str.startswith("B.1.1.529.2.86.1.1.11.1.1.1.3.8.1.2.1.2"),"Nextstrain_clade"] = "25F"  # NW.1.2
df.loc[df["unaliased"].str.startswith("XFC"),"Nextstrain_clade"] = "25G"  # XFC (recombinant)
df.loc[df["unaliased"].str.startswith("XFJ"),"Nextstrain_clade"] = "25H"  # XFJ (recombinant)
df.loc[df["unaliased"].str.startswith("B.1.1.529.3.2"),"Nextstrain_clade"] = "25I"  # BA.3.2

#%%
clade_muts = (
    df.groupby("Nextstrain_clade")
    .substitutions.apply(aggregate_mutations)
    .dropna()
    .astype(int)
)
#%%
clade_muts.rename("mut_count", inplace=True)
clade_muts
#%%
clade_count = df.groupby("Nextstrain_clade").count()
clade_count.rename(columns={"substitutions": "clade_count"}, inplace=True)
clade_count
# %%
#%%
mutations = pd.DataFrame(clade_muts).reset_index()
mutations.rename(columns={"level_1": "mutation"}, inplace=True)
mutations = mutations.join(clade_count, on="Nextstrain_clade")
mutations
#%%
mutations["proportion"] = mutations["mut_count"] / mutations["clade_count"]
mutations.sort_values(
    by=["Nextstrain_clade", "proportion"], ascending=False, inplace=True
)
mutations["genotype"] = mutations["mutation"].apply(lambda x: x[1:])
mutations["short_clade"] = mutations["Nextstrain_clade"].apply(lambda x: x[:3])
mutations

#%%
min_proportion = mutations["short_clade"].apply(lambda x: HIGH_THRESHOLD_PROPORTION if x in clades_with_high_proportion_threshold else MIN_PROPORTION)
# %%
# Choose which mutations to keep
relevant = mutations[
    (mutations["proportion"] > min_proportion)
]
relevant

# %% newly relevant
newly_relevant = relevant[mutations["proportion"] < 0.7]
print(newly_relevant[['Nextstrain_clade', 'mutation', 'proportion', 'mut_count']].to_csv())
# %%
mut_dict = {}
for mutation, row in relevant.groupby("genotype"):
    mut_dict[mutation] = (
        row[["short_clade", "mut_count"]]
        .sort_values(by="mut_count", ascending=False)["short_clade"]
        .to_list()
    )
mut_dict = dict(sorted(mut_dict.items(), key=lambda item: int(item[0][:-1])))
mut_dict

#%%
# Drop into /Users/corneliusromer/code/nextclade_data_workflows/sars-cov-2/defaults/labeled_muts.json
MUTS_TO_DROP = ["21T", "44T"]
for mut in MUTS_TO_DROP:
    mut_dict.pop(mut, None)
with open("defaults/labeled_muts.json", "w") as f_out:
    json.dump({"nucMutLabelMap": mut_dict}, f_out, indent=2)
#%%
def reverse_a_dict(dict_of_lists):
    result = {}
    for k, v in dict_of_lists.items():
        for x in v:
            result.setdefault(x, []).append(k)
    return result


with open("virus_properties.json", "w") as f_out:

    virus_json = {
        "schemaVersion": "1.10.0",
        "nucMutLabelMap": mut_dict,
        "nucMutLabelMapReverse": dict(sorted(reverse_a_dict(mut_dict).items())),
    }

    json.dump(virus_json, f_out, indent=2)

# %%
