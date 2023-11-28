#%%
from Bio import SeqIO
#%%
wuhan = str(SeqIO.read("/Users/corneliusromer/code/nextclade_data_workflows/sars-cov-2/profiles/clades/wuhan/reference.fasta", "fasta").seq)
xbb = str(SeqIO.read("/Users/corneliusromer/code/nextclade_data_workflows/sars-cov-2/profiles/clades/BA.2.86/reference.fasta", "fasta").seq)
ba2 = str(SeqIO.read("/Users/corneliusromer/code/nextclade_data_workflows/sars-cov-2/profiles/clades/21L/reference.fasta", "fasta").seq)
#%%
for i in range(len(wuhan)):
    if wuhan[i] != xbb[i] and wuhan[i] == ba2[i]:
        """
        "10029C": [
            "rev"
        ],
        """
        print(f'"{i+1}{wuhan[i]}": ["rev"],')
# %%
