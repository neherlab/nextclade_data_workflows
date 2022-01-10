# Produce files helpful for updating Nextclade's qc.json

rule select_frameshifts:
    input:
        "data/metadata_raw.tsv",
    output:
        "pre-processed/frameshifts.tsv"
    shell:
        """
        tsv-select -H -f frame_shifts {input} | tsv-summarize -H -w --group-by frame_shifts --count | keep-header -- sort -k2 -rn -t$'\\t' > {output} 
        """

rule select_stops:
    input:
        "data/metadata_raw.tsv",
    output:
        "pre-processed/aa_substitutions.tsv"
    shell:
        """
        tsv-select -H -f aaSubstitutions data/metadata_raw.tsv | grep '*' > {output} 
        """

rule filter_stops:
    input:
        "pre-processed/aa_substitutions.tsv"
    output:
        "pre-processed/stops_long.tsv"
    shell:
        """
        python scripts/filter_stops.py --input-file {input} --output-file {output}
        """

rule rank_stops:
    input:
        "pre-processed/stops_long.tsv",
    output:
        "pre-processed/stops.tsv"
    shell:
        """
        tsv-select -H -f stops {input} | tsv-summarize -H -w --group-by stops --count | keep-header -- sort -k2 -rn -t$'\\t' > {output} 
        """
