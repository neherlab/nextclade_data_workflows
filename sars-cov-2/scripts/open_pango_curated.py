#%%
from numpy import float64, int64
import pandas as pd
import datetime as dt
#%%
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_rows', 500 )
pd.set_option('display.width', 1000)

#%%
dtype_dict = {
    'strain': str,
    'virus': str,
    'gisaid_epi_isl': str,
    'genbank_accession': str,
    'sra_accession': str,
    'date': str,
    'region': str,
    'country': str,
    'division': str,
    'location': str,
    'region_exposure': str,
    'country_exposure': str,
    'division_exposure': str,
    'segment': str,
    'length': str,
    'host': str,
    'age': str,
    'sex': str,
    'Nextstrain_clade': str,
    'pango_lineage': str,
    'GISAID_clade': str,
    'originating_lab': str,
    'submitting_lab': str,
    'authors': str,
    'url': str,
    'title': str,
    'paper_url': str,
    'date_submitted': str,
    'sampling_strategy': str,
    'missing_data': str,
    'divergence': str,
    'nonACGTN': str,
    'rare_mutations': str,
    'snp_clusters': str,
    'QC_missing_data': str,
    'QC_mixed_sites': str,
    'QC_rare_mutations': str,
    'QC_snp_clusters': str,
    'clock_deviation': str
}
#%%
#download the master list of curated pangolin designations
#https://www.pangolin.org/downloads/curated_designations.txt
df_pango = pd.read_csv('https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.csv')
df_pango

#%%
# Load open data metadata into dataframe
df_open = pd.read_csv('data/metadata.tsv', sep='\t', dtype=dtype_dict)

#%%
df_gisaid = pd.read_csv('data/gisaid_metadata.tsv', sep='\t', dtype=dtype_dict)
#%%
#check for each curated strain whether it's found in the open metadata
#maybe need to speed up this process by using a dictionary of open metadata
df = df_pango.join(df_gisaid.set_index('strain'), on='taxon',lsuffix='', rsuffix='_gisaid',how='left')
df = df.join(df_open.set_index('strain'), on='taxon',lsuffix='', rsuffix='_open',how='left')
df['date_submitted'] = pd.to_datetime(df.date_submitted)
df
#%%
df.columns.tolist()
# %%
df = df.assign(in_open=df.date_submitted_open.apply(lambda x: 'Yes' if not type(x) == str else 'No'))
df = df.assign(in_gisaid=df.date_submitted.apply('notna'))
df.in_open.value_counts()
#%%
df.in_gisaid.value_counts()
#%%
# Output pivot table for each lineage with two columns one those sequences that are in open data and one for those that are not
df_lin = df.pivot_table(index='lineage', columns="in_open", aggfunc='count').taxon
#%%
df_lin.to_csv('temp')
# %%
df_lin['No'] = df_lin.No.fillna(0).astype(int)
df_lin['Yes'] = df_lin.Yes.fillna(0).astype(int)
df_lin['total'] = df_lin.sum(axis=1).fillna(0).astype(int)
df_lin['open_ratio'] = (df_lin['Yes'] / df_lin.total * 100).fillna(0).astype(int)
df_lin.to_csv('temp')
df_lin
#%%
# sort data frame by number of Yes downwards
df_lin.sort_values('Yes', ascending=False, inplace=True)
df_lin.to_csv('temp')
df_lin
# %%
# Add proportion of pango designation that is in the past 6 months
def recency(x):
    return 'NA' if not x['in_gisaid'] else('Recent' if x['date_submitted'] > dt.datetime(2021,4,1) else 'Old')

df['recency'] = df.apply(recency,axis=1)
#%%
df.recency.value_counts()
#%%
df_new = df.pivot_table(index='lineage', columns=["recency",'in_open'], values='taxon', aggfunc='count').fillna(0).astype(int)
df_new.sort_values(("Recent","Yes"), ascending=False, inplace=True)
df_new.to_csv('other_data/recent_curated_pango.tsv', sep='\t')
df_new
# %%
df