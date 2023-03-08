# Produce files helpful for updating Nextclade's qc.json


rule prepare_qc:
    input:
        "pre-processed/frameshifts.tsv",
        "pre-processed/stops.tsv",
        "pre-processed/stops_long.tsv",
        "pre-processed/stops.txt",
        "pre-processed/frameshifts.txt",


rule select_frameshifts:
    input:
        rules.download_metadata.output.metadata,
    output:
        "pre-processed/frameshifts_raw.tsv",
    shell:
        """
        zstdcat {input} | \
        tsv-select -H -f frame_shifts | tsv-filter -H --not-empty frame_shifts | tsv-summarize -H -w -x --group-by frame_shifts --count | keep-header -- sort -k2 -rn -t$'\\t' > {output} 
        """


rule exclude_impossible_frameshifts:
    input:
        frameshifts="pre-processed/frameshifts_raw.tsv",
        exclude="defaults/frameshifts_exclude.txt",
    output:
        "pre-processed/frameshifts.tsv",
    shell:
        """
        # Exclude all lines from the input file that are present in the exclude file
        grep -v -f {input.exclude} {input.frameshifts} > {output}
        """


rule select_stops:
    input:
        rules.download_metadata.output.metadata,
    output:
        "pre-processed/aa_substitutions.tsv",
    shell:
        """
        zstdcat {input} | \
        tsv-select -H -f aaSubstitutions | grep '*' > {output} 
        """


rule filter_stops:
    input:
        aa="pre-processed/aa_substitutions.tsv",
        exclude="defaults/stops_exclude.txt",
    output:
        "pre-processed/stops_long.tsv",
    shell:
        """
        python scripts/filter_stops.py --input-file {input.aa} --output-file /dev/stdout | \
        grep -v -f {input.exclude} > {output}
        """


rule rank_stops:
    input:
        "pre-processed/stops_long.tsv",
    output:
        "pre-processed/stops.tsv",
    shell:
        """
        tsv-select -H -f stops {input} | tsv-summarize -H -w --group-by stops --count | keep-header -- sort -k2 -rn -t$'\\t' > {output} 
        """


rule format_stops:
    input:
        "pre-processed/stops.tsv",
    output:
        "pre-processed/stops.txt",
    shell:
        """
        python scripts/common_stops.py --number 50 --input-file {input} >{output}
        """


rule format_frameshifts:
    input:
        "pre-processed/frameshifts.tsv",
    output:
        "pre-processed/frameshifts.txt",
    shell:
        """
        python scripts/common_frameshifts.py --number 20 --input-file {input} >{output}
        """
