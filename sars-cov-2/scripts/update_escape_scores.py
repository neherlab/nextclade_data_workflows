#%%
import json
from collections import defaultdict
import copy

#%%
# Load old virus_properties.json
old = json.load(open("sars-cov-2/old_virus_properties.json", "r"))
escape = json.load(open("sars-cov-2/escape_scores.json", "r"))
escape

#%%
old
#%%
new = copy.deepcopy(old)
new.keys()

#%%
new["phenotypeData"][0]["data"]
#%%
new["phenotypeData"][0]["data"] = escape["data"]
new["phenotypeData"][0]["data"]
#%%
with open("sars-cov-2/updated_virus_properties.json", "w") as f_out:

    json.dump(new, f_out, indent=2, sort_keys=False)

# %%
