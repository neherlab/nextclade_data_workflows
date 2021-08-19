#%%
import pandas as pd
import yaml

#%%
roots = ["CY115151", "CY114381", "CY121680", "CY115183"]

#%%
strainNames = []
for strain in ["A","B"]:
    df = pd.read_csv(
            f"pre-processed/metadata_{strain}_4.tsv",
            header=0,
            sep="\t",
        )
    strainNames.extend(df[df.ncbiAcc.isin(roots)].strainName.values)
strainNames
# %%

root = {}
for strain in ["A","B"]:
    for segment in range(1,8):
        df = pd.read_csv(
            f"pre-processed/metadata_{strain}_{segment}.tsv",
            header=0,
            sep="\t",
        )
        root[strain] = {segment : (df[df.strainName.isin(strainNames)].ncbiAcc.values)}
root
# %%
def acc(strainName, segment):
    df = pd.read_csv(
        f"pre-processed/metadata_{strainName[0]}_{segment}.tsv",
        header=0,
        sep="\t",
    )
    return (df[df.strainName ==strainName].ncbiAcc.values[0])
# %%
root = {}
for strain in strainNames:
    temp = {}
    for segment in range(1,9):
        temp[segment] = acc(strain,segment)
    root[strain] = temp
root
# %%
yaml.dump(root, open("acc.yaml", "w"))
# %%
