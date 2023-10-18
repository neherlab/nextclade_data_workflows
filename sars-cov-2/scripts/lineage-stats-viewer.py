# %%
import numpy as np
import pandas as pd

# %%
LINEAGE = "XCT"
# %%
rawdata = np.load("pre-processed/pango_matrix.npz")[LINEAGE]
# %%
array = rawdata.reshape(-1, 6)
# %%
df = pd.DataFrame(array, columns=["A", "C", "G", "T", "del", "N"])
df["pos"] = df.index + 1
# Make this column the index
df.set_index("pos", inplace=True)
df["argmax"] = df.idxmax(axis=1)
df["max"] = df.max(axis=1)

# %%
# Include column with each of reference and ba.2
# Load them with BioPython SeqIO from profiles/clades/wuhan/reference.fasta
import Bio.SeqIO

ref = Bio.SeqIO.read("profiles/clades/wuhan/reference.fasta", "fasta")
# Generate an array and append as new column
ref_array = np.array(list(str(ref.seq)))
df["ref"] = ref_array
df

# %%

df.to_csv(f"pre-processed/{LINEAGE}.csv", index=True)

# %%
# Then can look at it in visidata
