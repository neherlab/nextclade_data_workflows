#%%
import argparse, os, glob
from augur.io import open_file
from Bio import SeqIO, SeqFeature, Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

print(glob.glob("*"))
#%%
vic = pd.read_csv('pre-processed/vic.mutation_summary.tsv',sep='\t')
yam = pd.read_csv('pre-processed/yam.mutation_summary.tsv',sep='\t')
# %%
vic = vic.assign(vic_dist = vic.nucleotide.apply(lambda x: len(x.split(","))))
vic = vic.assign(yam_dist = yam.nucleotide.apply(lambda x: len(x.split(","))))
# yam.nucleotide.apply(lambda x: len(x.split(","))).hist()
# %%
vic.iloc[:,0] == yam.iloc[:,0]
# %%
plt.scatter(x=vic["vic_dist"],y=vic["yam_dist"])
plt.plot([0,100,200],[0,100,200])
plt.xlabel('vic distance')
plt.ylabel('yam distance')
# %%
plt.hist(vic.vic_dist - vic.yam_dist,bins=100)
# plt.ylim(0,10)
plt.yscale('log')
plt.ylabel('# of occurences')
plt.xlabel('vic distance - yam distance')
# %%
