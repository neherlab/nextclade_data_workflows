# Produce files helpful for updating Nextclade's qc.json


rule prepare_qc:
    input:
        "pre-processed/new_stops.ndjson",
        "pre-processed/new_frameshifts.ndjson",


rule recent_sequences:
    input:
        rules.download_metadata.output.metadata,
    output:
        "pre-processed/recent_sequences.tsv",
    shell:
        """
        zstdcat {input} | \
        tsv-filter -H --str-in-fld date:2023 | \
        tsv-select -H -f strain,aaSubstitutions,frame_shifts > {output}
        """


rule select_frameshifts:
    input:
        "pre-processed/recent_sequences.tsv",
    output:
        "pre-processed/frameshifts.tsv",
    shell:
        """
        zstdcat {input} \
        | tsv-select -H -f frame_shifts \
        | tsv-filter -H --not-empty frame_shifts \
        | tsv-filter -H --not-regex frame_shifts:'(ORF1a|ORF1b|N|M|S|E)' \
        | sed 's/,/\\n/g' \
        | tsv-summarize -H -w -x --group-by frame_shifts --count \
        | tsv-filter -H --ge count:20 \
        | keep-header -- sort -k2 -rn -t$'\\t' > {output} 
        """


rule select_stops:
    input:
        "pre-processed/recent_sequences.tsv",
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
        tsv-select -H -f stops {input} \
        | tsv-summarize -H -w --group-by stops --count \
        | tsv-filter -H --not-regex stops:'^(ORF1a|ORF1b|N|M|S|E):.*' \
        | tsv-filter -H --ge count:100 \
        | keep-header -- sort -k2 -rn -t$'\\t' \
        > {output} 
        """


rule format_stops:
    input:
        "pre-processed/stops.tsv",
    output:
        "pre-processed/common_stops.ndjson",
    shell:
        """
        python scripts/common_stops.py --number 100 --input-file {input} \
        | jq -c '.' >{output}
        """


rule format_frameshifts:
    input:
        "pre-processed/frameshifts.tsv",
    output:
        "pre-processed/common_frameshifts.ndjson",
    shell:
        """
        python scripts/common_frameshifts.py --number 500 --input-file {input} \
        | jq -c '.' >{output}
        """


rule old:
    input:
        qc_json="../../nextclade_data/data/datasets/sars-cov-2-21L/references/BA.2/versions/2023-09-21T12:00:00Z/files/qc.json",
    output:
        stops="pre-processed/old_stops.ndjson",
        frameshifts="pre-processed/old_frameshifts.ndjson",
    shell:
        """
        jq -c '.stopCodons.ignoredStopCodons[]' {input.qc_json} > {output.stops}
        jq -c '.frameShifts.ignoredFrameShifts[]' {input.qc_json} > {output.frameshifts}
        """


rule new_stops:
    """
    Only output lines from the new stops file that are not present in the old stops file
    """
    input:
        old="pre-processed/old_stops.ndjson",
        new="pre-processed/common_stops.ndjson",
    output:
        "pre-processed/new_stops.ndjson",
    run:
        import json
        from operator import itemgetter

        with open(input.old, "r") as f1, open(input.new, "r") as f2, open(
            output[0], "w"
        ) as f_out:
            old_stops = set(line.strip() for line in f1)
            new_stops = set(line.strip() for line in f2)
            diff = new_stops - old_stops
            diff_dict = [json.loads(stop) for stop in diff]
            sorted_stops = sorted(diff_dict, key=itemgetter("cdsName", "codon"))
            for stop in sorted_stops:
                f_out.write(json.dumps(stop) + "\n")


rule new_frameshifts:
    """
    Only output lines from the new stops file that are not present in the old stops file
    """
    input:
        old="pre-processed/old_frameshifts.ndjson",
        new="pre-processed/common_frameshifts.ndjson",
    output:
        "pre-processed/new_frameshifts.ndjson",
    run:
        import json


        def sorting_key(stop):
            gene_name = stop["cdsName"]
            begin_codon = stop["codonRange"]["begin"]
            end_codon = stop["codonRange"]["end"]
            return (gene_name, begin_codon, end_codon)


        with open(input.old, "r") as f1, open(input.new, "r") as f2, open(
            output[0], "w"
        ) as f_out:
            old = set(line.strip() for line in f1)
            new = set(line.strip() for line in f2)
            print(old)
            print(new)
            diff = new - old
            diff_dict = [json.loads(item) for item in diff]
            sorted_items = sorted(diff_dict, key=sorting_key)
            for item in sorted_items:
                f_out.write(json.dumps(item) + "\n")
