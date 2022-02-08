#%%
import json
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

#%%
# Defines minimum proportions and counts for relevant mutations to keep
MIN_PROPORTION = 0.3
MIN_COUNT = 100000
#%%
# aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.gz .
df = pd.read_csv(
    # "metadata.tsv.gz",
    "s3://nextstrain-ncov-private/metadata.tsv.gz",
    sep="\t",
    usecols=["Nextstrain_clade", "substitutions"],
    parse_dates=False,
).dropna()
df
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
# %%
# Choose which mutations to keep
relevant = mutations[
    (mutations["proportion"] > MIN_PROPORTION) | (mutations["mut_count"] > MIN_COUNT)
]
relevant
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
