#%%
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_rows', 5) 
pd.set_option('display.width', 1000)
#%%
#download the master list of curated pangolin designations
#https://www.pangolin.org/downloads/curated_designations.txt
df = pd.read_csv('https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.metadata.csv')

#%%
# Load open data metadata into dataframe
df_open = pd.read_csv('data/metadata.tsv', sep='\t')
#%%
#check for each curated strain whether it's found in the open metadata
#maybe need to speed up this process by using a dictionary of open metadata
df = df.join(df_open.set_index('strain'), on='sequence_name',lsuffix='_pango', rsuffix='_open',how='left')

# %%
df = df.assign(in_open=df.date_submitted.notnull().apply(lambda x: 'Yes' if x else 'No'))
#%%
# Output pivot table for each lineage with two columns one those sequences that are in open data and one for those that are not
df_lin = df.pivot_table(index='lineage', columns="in_open", values='sequence_name', aggfunc='count')
df_lin
# %%
df_lin['total'] = df_lin.sum(axis=1)
df_lin['open_ratio'] = df_lin['Yes'] / df_lin.total
#%%
# sort data frame by number of Yes downwards
df_lin.sort_values('Yes', ascending=False, inplace=True)
df_lin
#%%
df_lin.sort_index
#%%
# output to tsv
df_lin.to_csv('other_data/curated_pango_open_stats.tsv', sep='\t')
# %%
df_lin
# %%
