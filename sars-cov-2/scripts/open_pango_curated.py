#%%
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
df = pd.read_csv('https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.csv')
df

#%%
# Load open data metadata into dataframe
df_open = pd.read_csv('data/metadata.tsv', sep='\t', dtype=dtype_dict)

#%%
df_gisaid = pd.read_csv('data/gisaid_metadata.tsv', sep='\t', dtype=dtype_dict)
#%%
#check for each curated strain whether it's found in the open metadata
#maybe need to speed up this process by using a dictionary of open metadata
df = df.join(df_gisaid.set_index('strain'), on='taxon',lsuffix='', rsuffix='_gisaid',how='left')
df = df.join(df_open.set_index('strain'), on='taxon',lsuffix='', rsuffix='_open',how='left')
df
#%%
df.columns.tolist()
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
df_gisaid[df_gisaid['strain'].str.contains('ONP180821NB32')]
#%%
# sort data frame by number of Yes downwards
df_lin.sort_values('Yes', ascending=False, inplace=True)
df_lin
#%%
# output to tsv
df_lin.to_csv('other_data/curated_pango_open_stats.tsv', sep='\t')
# %%
df_lin
# %%
# Add proportion of pango designation that is in the past 6 months
df['recency'] = df.sample_date.apply(lambda x: 'Recent' if x > dt.datetime(2021,4,1) else 'Old')
df
#%%
df_new = df.pivot_table(index='lineage', columns=["recency","in_open"], values='sample_date', aggfunc='count')
df_new.to_csv('other_data/recent_curated_pango.tsv', sep='\t')
# %%
for key, value in df.dtypes.items():
    print(key, value)
# %%
for row in df_new.iterrows():
    print(row)
# %%
