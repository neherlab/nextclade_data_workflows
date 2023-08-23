"""
This script makes updates to virus_properties when a new clade is added
New clades require new characteristic mutations
But we need to make sure we don't throw out other important information in virus properties
Hence this script loads both old and new virus_properties.json and makes the necessary updates
"""
#%%
import copy
import json

#%%
path_old = "/Users/corneliusromer/code/nextclade_data/data/datasets/sars-cov-2-21L/references/BA.2/versions/2023-08-09T12:00:00Z/files/virus_properties.json"
# path_old = "/Users/corneliusromer/code/nextclade_data/data/datasets/sars-cov-2/references/MN908947/versions/2023-08-09T12:00:00Z/files/virus_properties.json"
path_new = "/Users/corneliusromer/code/nextclade_data_workflows/sars-cov-2/virus_properties.json"
# path_out = "virus_properties_updated.json"
path_out = path_old
# Load old virus_properties.json
old = json.load(open(path_old, "r"))
new = json.load(open(path_new, "r"))
old_copy = copy.deepcopy(old)

# Need to keep the extra data
# Pull in "rev" data
# Replace those mutations that are now rev
#%%
#%%
for k,v in new["nucMutLabelMap"].items():
    old["nucMutLabelMap"][k] = copy.deepcopy(v) #list(filter(lambda x: x != "rec", v))
    print(k)

#%%
# Sort old["nucMutLabelMap"]
old["nucMutLabelMap"] = dict(sorted(old["nucMutLabelMap"].items(), key=lambda item: int(item[0][:-1])))
#%%
# Don't use these mutations as "characteristic" as they are flip floppy
for k in ["21T", "44T"]:
    try:
        old["nucMutLabelMap"].pop(k)
    except KeyError:
        pass

def reverse_a_dict(dict_of_lists):
    result = {}
    for k, v in dict_of_lists.items():
        for x in v:
            result.setdefault(x, []).append(k)
    return result

old["nucMutLabelMapReverse"] = dict(sorted(reverse_a_dict(old["nucMutLabelMap"]).items()))


with open(path_out, "w") as f_out:
    json.dump(old, f_out, indent=None, sort_keys=False)

# %%
