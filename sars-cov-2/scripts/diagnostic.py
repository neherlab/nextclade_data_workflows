"""
Use nextclade QC to produce a list of sequences to be excluded.
"""
import argparse
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def datestr_to_ordinal(x, minus_weeks=0):
    try:
        return (datetime.strptime(x,"%Y-%m-%d") - timedelta(weeks=minus_weeks)).toordinal()
    except:
        return np.nan

def earliest_clade_date(Nextstrain_clade, clade_emergence_dates_filename, window_weeks=2):
    clade_dates = pd.read_csv(clade_emergence_dates_filename, index_col="Nextstrain_clade", sep='\t')
    try:
        return datestr_to_ordinal(clade_dates.loc[Nextstrain_clade]['first_sequence'], minus_weeks=window_weeks)
    except:
        return np.nan

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", type=str, required=True, help="metadata tsv file")
    parser.add_argument("--clade_emergence_dates", type=str, default="defaults/clade_emergence_dates.tsv", help="tsv file with two columns: Nextstrain_clade name and first known sequence for that clade.")
    parser.add_argument("--clock-filter-lower-limit", type=float, default=-10, help="tsv file with two columns: Nextstrain_clade name and first known sequence for that clade.")
    parser.add_argument("--clock-filter-recent", type=float, default=20, help="max allowed clock deviation for recently submitted sequences")
    parser.add_argument("--clock-filter", type=float, default=15, help="max allowed clock deviation for non-recent sequences")
    parser.add_argument("--snp-clusters", type=int, default=1, help="max allowed SNP clusters (called by nextclade)")
    parser.add_argument("--rare-mutations", type=int, default=35, help="max allowed private mutations (called by nextclade)")
    parser.add_argument("--clock-plus-rare", type=int, default=40, help="maximal allowed clock deviation + rare mutations (called by nextclade)")
    parser.add_argument("--clade-emergence-window", type=int, default=2, help="number of weeks before official emergence of clade at which sequences can safely be excluded")
    parser.add_argument("--output-exclusion-list", type=str, required=True, help="Output to-be-reviewed addition to exclude.txt")
    parser.add_argument("--output-exclusion-reasons", type=str, help="Output reasons for exclusion as tsv")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t')

    check_recency = "date_submitted" in metadata.columns
    check_clade_dates = "Nextstrain_clade" in metadata.columns
    check_clock_deviation = "clock_deviation" in metadata.columns
    check_rare_mutations = "rare_mutations" in metadata.columns
    check_snp_clusters = "snp_clusters" in metadata.columns
    check_QC_mixed_sites= "QC_mixed_sites" in metadata.columns

    if check_recency:
        recency_cutoff = (datetime.today() - timedelta(weeks=6)).toordinal()
        recent_sequences = metadata.date_submitted.apply(lambda x: datestr_to_ordinal(x)>recency_cutoff)
    else:
        print("Skipping QC steps which rely on submission recency, as metadata is missing 'date_submitted'")

    if check_clade_dates:
        dates = metadata.date.apply(lambda x: datestr_to_ordinal(x))
        clade_dates = metadata.Nextstrain_clade.apply(lambda x: earliest_clade_date(x, args.clade_emergence_dates, window_weeks=args.clade_emergence_window))
    else:
        print("Skipping QC steps which rely on clade-date combinations, as metadata is missing 'Nextstrain_clade'")

    if check_clock_deviation:
        clock_deviation = np.array([float(x) if isfloat(x) else np.nan for x in metadata.clock_deviation])
    else:
        print("Skipping QC steps which rely on clock deviation, as metadata is missing 'clock_deviation'")

    if check_rare_mutations:
        rare_mutations = np.array([float(x) if isfloat(x) else np.nan for x in metadata.rare_mutations])
    else:
        print("Skipping QC steps which rely on rare mutations, as metadata is missing 'rare_mutations'")

    if check_snp_clusters:
        snp_clusters = np.array([float(x) if isfloat(x) else np.nan for x in metadata.snp_clusters])
    else:
        print("Skipping QC steps which rely on SNP clusters, as metadata is missing 'snp_clusters'")

    if check_QC_mixed_sites:
        QC_mixed_sites = metadata['QC_mixed_sites']
    else:
        print("Skipping QC steps which rely on mixed sites, as metadata is missing 'QC_mixed_sites'")

    rules = {
        "check_clade_dates": (dates<clade_dates) if check_clade_dates else None,
        "clock_filter_recent": np.abs(clock_deviation)>args.clock_filter_recent if check_clock_deviation else None,
        "clock_filter_old": (np.abs(clock_deviation)>args.clock_filter)&(~recent_sequences) if check_recency and check_clock_deviation else None,
        "clock_filter_lower_limit": clock_deviation<args.clock_filter_lower_limit if check_clock_deviation else None,
        "clock_plus_rare": np.abs(clock_deviation+rare_mutations)>args.clock_plus_rare if check_clock_deviation else None,
        "rare_mutation": rare_mutations>args.rare_mutations if check_rare_mutations else None,
        "snp_clusters": snp_clusters>args.snp_clusters if check_snp_clusters else None,
        "QC_mixed_sites": QC_mixed_sites == "bad" if check_QC_mixed_sites else None
    }

    exclude_meta = metadata.copy()
    exclude_meta['reason'] = ""

    to_exclude = np.zeros_like(clock_deviation, dtype=bool)
    for rule_name, rule_array in rules.items():
        if rule_array is not None:
            to_exclude = np.logical_or(to_exclude, rule_array)
            print(f"Excluded because of:   {rule_name:<20} {np.sum(rule_array):7}")
            exclude_meta.loc[rule_array,'reason'] += rule_name + ","
    exclude_meta['reason'] = exclude_meta['reason'].apply(lambda x: x.strip(","))

    if args.output_exclusion_reasons:
        with open(args.output_exclusion_reasons, 'w') as f:
            exclude_meta[to_exclude].to_csv(f, sep='\t', index=False, columns=["strain","reason","clock_deviation","date","date_submitted","country","Nextstrain_clade","pango_lineage","missing_data","divergence","rare_mutations","submitting_lab"])

    # write out file with sequences flagged for exclusion
    with open(args.output_exclusion_list, 'w') as excl:
        for s in metadata.loc[to_exclude,'strain']:
            excl.write(f'{s}\n')
