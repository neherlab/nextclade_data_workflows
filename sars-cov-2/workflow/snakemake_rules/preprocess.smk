'''
This part of the workflow downloads files from S3

  - "data/{origin}/sequences.fasta.xz"
  - "data/{origin}/metadata.tsv"

and produces

  - pre-processed/filtered.fasta.xz
  - pre-processed/metadata.tsv

'''

import os
localrules: download_sequences, download_metadata, download_exclude, download_clades, preprocess, download_color_ordering, download_curated_pango

rule preprocess:
    input:
        sequences = "data/sequences.fasta.xz",
        metadata = "data/metadata.tsv",
        sequence_index = "pre-processed/sequence_index.tsv",
        problematic_exclude = "pre-processed/problematic_exclude.txt",
    params:
        slack_hook = config.get('slackHook',"google.com")
    shell:
        """
        curl -X POST -H 'Content-type: application/json' \
        --data '{{"text":"Preprocessing done"}}' \
        {params.slack_hook}
        """

def _infer_decompression(input):
    """
    Returns a shell command to decompress the piped stream,
    which will itself produce a stream of decompressed data to stdout.
    If no decompression is needed, returns `cat`.
    NOTE: a lot of this will become unnecessary once `augur` handles
    compressed sequence inputs.
    """
    if input.endswith(".xz"):
        return "xz -dcq"
    if input.endswith(".gz"):
        return "gunzip -cq"
    return "cat"

rule download_sequences:
    message: "Downloading sequences from {params.address} -> {output[0]}"
    params:
        address = lambda w: config['origins']['sequences']
    output:
        "data/sequences.fasta.xz"
    shell: "curl {params.address} -o {output}"

rule download_metadata:
    message: "Downloading metadata from {params.address} -> {output}"
    params:
        deflate = lambda w: _infer_decompression(config['origins']['metadata']),
        address = lambda w: config['origins']['metadata']
    output:
        metadata = "data/metadata_raw.tsv"
    shell: "curl {params.address} | {params.deflate} > {output:q}"

rule fix_metadata:
    input: "data/metadata_raw.tsv"
    output: "data/metadata.tsv"
    shell: 
        """
        awk  -F'\t' 'BEGIN {{OFS = FS}} {{if ($NF=="?") $NF="-inf"; print}}' {input} >{output}
        """

rule download_exclude:
    message: "Downloading exclude from {params.source} -> {output}"
    output:
        "data/exclude.txt"
    params:
        source = lambda w: config["origins"]['exclude']
    shell: "curl {params.source} -o {output}"

rule download_clades:
    message: "Downloading clade definitions from {params.source} -> {output}"
    output:
        "builds/clades.tsv"
    params:
        source = config["data_source"]["clades"]
    shell: "curl {params.source} -o {output}"

rule download_color_ordering:
    message: "Downloading clade definitions from {params.source} -> {output}"
    output:
        "builds/color_ordering.tsv"
    params:
        source = config["data_source"]["color_ordering"]
    shell: "curl {params.source} -o {output}"

rule download_lat_longs:
    message: "Downloading clade definitions from {params.source} -> {output}"
    output:
        "builds/lat_longs.tsv"
    params:
        source = config["data_source"]["lat_longs"]
    shell: "curl {params.source} -o {output}"

rule download_curated_pango:
    output:
        "pre-processed/pango_raw.csv"
    params:
        source = config["data_source"]["pango"]
    shell: "curl {params.source} -o {output}"

rule strip_pango_strain_names:
    input:
        "pre-processed/pango_raw.csv"
    output:
        "pre-processed/pango.csv"
    shell: "awk -F',' '{{print $1}}' {input} | tail -n+2 >{output}"

rule diagnostic:
    message: "Scanning metadata {input.metadata} for problematic sequences. Removing sequences with >{params.clock_filter} deviation from the clock and with more than {params.snp_clusters}."
    input:
        metadata = "data/metadata.tsv"
    output:
        to_exclude = "pre-processed/problematic_exclude.txt"
    params:
        clock_filter_floor = -8,
        clock_filter_ceil = 20,
        snp_clusters = 1,
        rare_mutations = 30,
        clock_plus_rare = 100,
    log:
        "logs/diagnostics.txt"
    benchmark:
        "benchmarks/diagnostics.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000
    shell:
        """
        python3 scripts/diagnostic.py \
            --metadata {input.metadata} \
            --clock-filter-floor {params.clock_filter_floor} \
            --clock-filter-ceil {params.clock_filter_ceil} \
            --rare-mutations {params.rare_mutations} \
            --clock-plus-rare {params.clock_plus_rare} \
            --snp-clusters {params.snp_clusters} \
            --output-exclusion-list {output.to_exclude} 2>&1 | tee {log}
        """

rule index_sequences:
    message:
        """
        Index sequence composition for faster filtering.
        """
    input:
        sequences = "data/sequences.fasta.xz"
    output:
        sequence_index = "pre-processed/sequence_index.tsv"
    log:
        "logs/index_sequences.txt"
    benchmark:
        "benchmarks/index_sequences.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index} 2>&1 | tee {log}
        """

rule open_pango:
    # include only sequences that are in pango.csv using augur filter
    input:
        sequences = "data/sequences.fasta.xz",
        pango = "pre-processed/pango.csv",
        sequence_index = "pre-processed/sequence_index.tsv",
        metadata = "data/metadata.tsv",
    output: 
        sequences = "pre-processed/open_pango.fasta.xz",
        metadata = "pre-processed/open_pango_metadata.tsv",
        strains = "pre-processed/open_pango_strains.txt",
    log:
        "logs/open_pango.txt"
    benchmark:
        "benchmarks/open_pango.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --exclude-all \
            --include {input.pango} \
            --output-metadata {output.metadata} \
            --output-strains {output.strains} \
            --output-sequences {output.sequences} \
            2>&1 | tee {log} 
        """