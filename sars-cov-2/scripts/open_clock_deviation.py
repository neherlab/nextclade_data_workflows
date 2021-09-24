#%%
import pandas as pd
import datetime as dt
import numpy as np
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
# Load open data metadata into dataframe
df = pd.read_csv('data/gisaid_metadata.tsv', sep='\t', dtype=dtype_dict)
#%% 
df['cd'] = pd.to_numeric(df['clock_deviation'], errors='coerce')
#%%
#%%
df.dropna(subset=['cd'], inplace=True)
#%%
df['clock_deviation_bin'] = pd.cut(df['clock_deviation'], bins=range(20,100,1))
#%%
#%%
# Make histogram of clock deviation with a bin width of 1, for values between -20 and 20
def clock_deviation_histogram(df):
    df['clock_deviation'] = pd.to_numeric(df['clock_deviation'], errors='coerce')
    df['clock_deviation_bin'] = pd.cut(df['clock_deviation'], bins=range(20,100,1))
    return df.groupby(['clock_deviation_bin'])['clock_deviation'].count()
#%%
def quantile_q(q):
    def quantile(x):
        return round(np.percentile(x, q))
    quantile.__name__ = 'q' + str(q)
    return quantile

def interq(a,b):
    def interq(x):
        return round(-np.percentile(x, a) + np.percentile(x, b))
    interq.__name__ = 'iq' + f"{a}-{b}"
    return interq
#%%
def clock_pivot(df,index):
    df_temp = df.copy()
    # df_temp[index] = df[index].apply(lambda x: x[:50] if (type(x) == str) else x)
    df_pivot = df_temp.pivot_table(index= [index],values=['cd'],aggfunc=['count','min',quantile_q(5),quantile_q(25),quantile_q(50),quantile_q(75),quantile_q(95),'max'])
    # df_pivot.sort_values(by=[('count','clock_deviation')], ascending=False, inplace=True)
    # df_pivot.to_csv('other_data/clock_by_author.tsv', sep='\t')
    dfc = df_pivot[df_pivot[('count','cd')] > 100].sort_values(by=[('q25','cd')], ascending=True)
    dfc.to_csv(f"other_data/cd_by_{index}.tsv", sep='\t')
#%%