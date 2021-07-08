#%%
import treetime 
import treetime.utils
import Bio.Phylo
import Bio.Phylo.Newick as nwk
from typing import List
import numpy as np
import pandas as pd
import json
#%%
dates = treetime.utils.parse_dates("pre-processed/metadata_enriched.tsv")
myTree = treetime.TreeTime(dates=dates,tree="build/tree.nwk",gtr="JC69",seq_len=1600)
res = myTree.clock_filter(plot=True,n_iqd=1.5)
n_bad_after = [n.name for n in myTree.tree.get_terminals() if n.bad_branch]
# %%
res
# %%
tree : nwk.Tree = myTree.tree  
# %%
internal = [clade.branch_length for clade in tree.get_nonterminals()]
# %%
import matplotlib.pyplot as plt
plt.hist(internal)
# %%
clades : List[nwk.Clade] = tree.get_nonterminals()
# %%
# plt.scatter([clade.branch_length for clade in clades],[1/clade.count_terminals() for clade in clades])
# plt.hist([np.log((clade.branch_length+0.001)/clade.count_terminals()) for clade in clades])
plt.hist([np.log10((clade.branch_length+0.1)/clade.count_terminals()) for clade in clades],bins=50)
# %%
# %%

t = tree.get_terminals()
t
#%%
terminals = pd.DataFrame({"strain": [clade.name for clade in t], "clade": t})

# %%
terminals['length'] = terminals.clade.apply(lambda x: (x.branch_length+0.00))
# %%
terminals.sort_values(by="length")[-30:]
# %%
terminals['bad'] = False

# %%
# %%
terminals['bad'] = terminals.strain.isin(n_bad_after)

# %%
terminals[terminals.bad]
#%%
# %%
for node in n_bad_after:
    tree.prune(node)

# %%
Bio.Phylo.write(tree, "build/pruned_tree.nwk", 'newick')
# %%
dev = []
node_data = {'nodes':{}}
for name,val in res.items():
    node_data['nodes'][name.name] = {"Clock deviation": val}
    dev.append(val)

# %%
with open("build/clock_deviation.json",'w') as fp:
    json.dump(node_data,fp)
# %%
node_data = {'nodes':{}}
for idx, row in terminals.iterrows():
    node_data['nodes'][row.strain] = {"Terminal Length": row.length}
with open("build/terminal_branch_length.json",'w') as fp:
    json.dump(node_data,fp)
# %%
# %%
