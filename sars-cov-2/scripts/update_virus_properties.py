"""
This script makes updates to virus_properties when a new clade is added
New clades require new characteristic mutations
But we need to make sure we don't throw out other important information in virus properties
Hence this script loads both old and new virus_properties.json and makes the necessary updates
"""
#%%
import copy
import json

path_old = "/Users/corneliusromer/code/nextclade_data/data/nextstrain/sars-cov-2/wuhan-hu-1/orfs/pathogen.json"
# path_old = "/Users/corneliusromer/code/nextclade_data/data/datasets/sars-cov-2/references/MN908947/versions/2023-08-09T12:00:00Z/files/virus_properties.json"
path_new = "/Users/corneliusromer/code/nextclade_data_workflows/sars-cov-2/virus_properties.json"

def update_pathogen_json(path_old, path_new):
    # path_out = "virus_properties_updated.json"
    path_out = path_old
    # Load old virus_properties.json
    old_full = json.load(open(path_old, "r"))
    old = old_full["mutLabels"]
    new = json.load(open(path_new, "r"))
    old_copy = copy.deepcopy(old)

    # Need to keep the extra data
    # Pull in "rev" data
    # Replace those mutations that are now rev
    
    
    # Leave reversions in place - just add new mutations
    for k,v in new["nucMutLabelMap"].items():
        old["nucMutLabelMap"][k] = copy.deepcopy(v) #list(filter(lambda x: x != "rec", v))
        print(k)

    
    # Sort old["nucMutLabelMap"]
    old["nucMutLabelMap"] = dict(sorted(old["nucMutLabelMap"].items(), key=lambda item: int(item[0][:-1])))
    
    # Don't use these mutations as "characteristic" as they are flip floppy
    for k in ["21T", "44T"]:
        try:
            old["nucMutLabelMap"].pop(k)
        except KeyError:
            pass

    # def reverse_a_dict(dict_of_lists):
    #     result = {}
    #     for k, v in dict_of_lists.items():
    #         for x in v:
    #             result.setdefault(x, []).append(k)
    #     return result

    # No longer needed
    # old["nucMutLabelMapReverse"] = dict(sorted(reverse_a_dict(old["nucMutLabelMap"]).items()))

    old_full["mutLabels"] = old

    with open(path_out, "w") as f_out:
        json.dump(old_full, f_out, indent=2, sort_keys=False)

#%%    
path_old = "/Users/corneliusromer/code/nextclade_data/data/nextstrain/sars-cov-2/wuhan-hu-1/orfs/pathogen.json"
path_new = "/Users/corneliusromer/code/nextclade_data_workflows/sars-cov-2/virus_properties.json"
update_pathogen_json(path_old, path_new)

#
paths = [f"/Users/corneliusromer/code/nextclade_data/data/nextstrain/sars-cov-2/{d}/pathogen.json" for d in ["wuhan-hu-1/orfs", "wuhan-hu-1/proteins", "BA.2", "BA.2.86", "XBB"]]
for path in paths:
    update_pathogen_json(path, path_new)

# %%
