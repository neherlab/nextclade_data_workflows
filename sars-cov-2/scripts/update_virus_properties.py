#%%
import json
from collections import defaultdict
import copy

#%%
# Load old virus_properties.json
old = json.load(open("sars-cov-2/old_virus_properties.json", "r"))
new = json.load(open("sars-cov-2/virus_properties.json", "r"))

#%%
for k,v in new["nucMutLabelMap"].items():
    old["nucMutLabelMap"][k] = list(filter(lambda x: x != "rec", v))

for k in ["21T", "44T"]:
    old["nucMutLabelMap"].pop(k)

def reverse_a_dict(dict_of_lists):
    result = {}
    for k, v in dict_of_lists.items():
        for x in v:
            result.setdefault(x, []).append(k)
    return result

old["nucMutLabelMapReverse"] = dict(sorted(reverse_a_dict(old["nucMutLabelMap"]).items()))


with open("sars-cov-2/updated_virus_properties.json", "w") as f_out:

    json.dump(old, f_out, indent=None, sort_keys=False)

# %%
