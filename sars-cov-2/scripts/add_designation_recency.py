"""
Reads in designation_dates.tsv which contains `lineage,designation_date` pairs
Adds a new column that calculates recency in categories:
- 0-2 days
- 3-7 days
- 1-2 weeks
- 2-3 weeks
- 3-4 weeks
- 4-8 weeks
- >8 weeks

Usage:
python3 add_designation_recency.py <path/to/lineages.tsv> <path/to/lineages_with_recency.tsv>

Using typer to create a CLI
"""
import typer
import pandas as pd

def add_designation_recency(lineages_tsv: str, output_tsv: str):
    df = pd.read_csv(lineages_tsv, sep='\t')
    df['recency'] = (pd.to_datetime('today') - pd.to_datetime(df['designation_date'])).dt.days
    df.loc[df['recency'] <= 2, 'recency_category'] = '0-2 days'
    df.loc[(df['recency'] > 2) & (df['recency'] <= 7), 'recency_category'] = '3-7 days'
    df.loc[(df['recency'] > 7) & (df['recency'] <= 14), 'recency_category'] = '1-2 weeks'
    df.loc[(df['recency'] > 14) & (df['recency'] <= 21), 'recency_category'] = '2-3 weeks'
    df.loc[(df['recency'] > 21) & (df['recency'] <= 28), 'recency_category'] = '3-4 weeks'
    df.loc[(df['recency'] > 28) & (df['recency'] <= 56), 'recency_category'] = '4-8 weeks'
    df.loc[df['recency'] > 56, 'recency_category'] = '>8 weeks'
    df = df.drop(columns=['recency'])
    df.rename(columns={'recency_category': 'designation_recency'}, inplace=True)
    df.to_csv(output_tsv, sep='\t', index=False)

if __name__ == '__main__':
    typer.run(add_designation_recency)
