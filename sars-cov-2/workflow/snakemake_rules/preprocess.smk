'''
This part of the workflow downloads files from S3

  - "data/{origin}/sequences.fasta.xz"
  - "data/{origin}/metadata.tsv"

and produces

  - pre-processed/filtered.fasta.xz
  - pre-processed/metadata.tsv

'''

import os
localrules: download_sequences, download_metadata, download_exclude, download_clades, preprocess

rule preprocess:
    input:
        sequences = "pre-processed/filtered.fasta.xz",
        metadata = "pre-processed/metadata.tsv",
        sequence_index = "pre-processed/sequence_index.tsv"
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
        address = lambda w: config['origins'][w.origin]['sequences']
    output:
        "data/{origin}/sequences.fasta.xz"
    shell: "conda activate nextstrain; aws s3 cp {params.address} {output}"

rule download_metadata:
    message: "Downloading metadata from {params.address} -> {output}"
    params:
        deflate = lambda w: _infer_decompression(config['origins'][w.origin]['metadata']),
        address = lambda w: config['origins'][w.origin]['metadata']
    output:
        metadata = "data/{origin}/metadata.tsv"
    shell: "conda activate nextstrain; aws s3 cp {params.address} - | {params.deflate} {input} > {output:q}"

rule download_exclude:
    message: "Downloading exclude from {params.source} -> {output}"
    output:
        "data/{origin}/exclude.txt"
    params:
        source = lambda w: config["origins"][w.origin]['exclude']
    shell: "curl {params.source} -o {output}"

rule download_clades:
    message: "Downloading clade definitions from {params.source} -> {output}"
    output:
        "data/clades.tsv"
    params:
        source = config["data_source"]["clades"]
    shell: "curl {params.source} -o {output}"

rule prealign:
    message:
        """
        Aligning sequences to {input.reference}
            - gaps relative to reference are considered real
        """
    input:
        sequences = "data/{origin}/sequences.fasta.xz",
        genemap = config["files"]["annotation"],
        reference = config["files"]["alignment_reference"]
    output:
        alignment = "pre-processed/{origin}/alignment.fasta.xz",
        insertions = "pre-processed/{origin}/insertions.tsv",
        translations = expand("pre-processed/{{origin}}/translations/seqs.gene.{gene}.fasta.xz", gene=config.get('genes', ['S']))
    params:
        outdir = "pre-processed/{origin}/translations",
        genes = ','.join(config.get('genes', ['S'])),
        basename = "seqs",
        tmp_alignment = "pre-processed/{origin}/alignment.fasta",
        deflate = lambda w: _infer_decompression(".xz")
    log:
        "logs/prealign_{origin}.txt"
    benchmark:
        "benchmarks/prealign_{origin}.txt"
    conda: config["conda_environment"]
    threads: 8
    resources:
        mem_mb = 3000,
    shell:
        """
	conda activate nextstrain;
        {params.deflate} {input.sequences} | nextalign \
            --jobs={threads} \
            --reference {input.reference} \
            --genemap {input.genemap} \
            --genes {params.genes} \
            --sequences /dev/stdin \
            --output-dir {params.outdir} \
            --output-basename {params.basename} \
            --output-fasta {params.tmp_alignment} \
            --output-insertions {output.insertions} &&\
	    xz -2 {params.tmp_alignment} &&\
        xz -2 {params.outdir}/*fasta
        """

rule diagnostic:
    message: "Scanning metadata {input.metadata} for problematic sequences. Removing sequences with >{params.clock_filter} deviation from the clock and with more than {params.snp_clusters}."
    input:
        metadata = "data/{origin}/metadata.tsv"
    output:
        to_exclude = "pre-processed/{origin}/problematic_exclude.txt"
    params:
        clock_filter = 20,
        snp_clusters = 1,
        rare_mutations = 100,
        clock_plus_rare = 100,
    log:
        "logs/diagnostics_{origin}.txt"
    benchmark:
        "benchmarks/diagnostics_{origin}.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000
    shell:
        """
	conda activate nextstrain;
        python3 scripts/diagnostic.py \
            --metadata {input.metadata} \
            --clock-filter {params.clock_filter} \
            --rare-mutations {params.rare_mutations} \
            --clock-plus-rare {params.clock_plus_rare} \
            --snp-clusters {params.snp_clusters} \
            --output-exclusion-list {output.to_exclude} 2>&1 | tee {log}
        """

rule filter:
    message:
        """
        Filtering alignment {input.sequences} -> {output.sequences}
          - excluding strains in {input.exclude} {input.problematic}
          - including strains in {input.include}
        """
    input:
        sequences = "pre-processed/{origin}/alignment.fasta.xz",
        metadata = "data/{origin}/metadata.tsv",
        include = "defaults/include.txt",
        exclude = "data/{origin}/exclude.txt",
        problematic = "pre-processed/{origin}/problematic_exclude.txt"
    output:
        sequences = "pre-processed/{origin}/filtered.fasta.xz"
    log:
        "logs/filtered{origin}.txt"
    benchmark:
        "benchmarks/filter{origin}.txt"
    params:
        filter_arguments = lambda w: config["origins"][w.origin].get("filters",""),
        tmp_alignment = "pre-processed/{origin}/filtered.fasta"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
	conda activate nextstrain;
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --exclude {input.exclude} {input.problematic} \
            {params.filter_arguments} \
            --output {params.tmp_alignment} 2>&1 | tee {log};
        xz -2 {params.tmp_alignment}
        """


rule combine_bulk_sequences:
    input:
        sequences = [f"pre-processed/{origin}/filtered.fasta.xz" for origin in config["origins"]],
        mutation_summary = [f"pre-processed/{origin}/mutation_summary.tsv" for origin in config["origins"]]
    output:
        rules.preprocess.input.sequences
    run:
        if len(input.sequences)==1:
            shell(f"cp {input.sequences} {output}")

rule combine_bulk_metadata:
    input:
        [f"data/{origin}/metadata.tsv" for origin in config["origins"]]
    output:
        rules.preprocess.input.metadata
    run:
        if len(input)==1:
            shell(f"cp {input} {output}")

rule index_sequences:
    message:
        """
        Index sequence composition for faster filtering.
        """
    input:
        sequences = rules.combine_bulk_sequences.output
    output:
        sequence_index = rules.preprocess.input.sequence_index
    log:
        "logs/index_sequences.txt"
    benchmark:
        "benchmarks/index_sequences.txt"
    conda: config["conda_environment"]
    shell:
        """
	conda activate nextstrain;
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index} 2>&1 | tee {log}
        """


rule mutation_summary:
    message: "Summarizing {input.alignment}"
    input:
        alignment = rules.prealign.output.alignment,
        insertions = rules.prealign.output.insertions,
        translations = rules.prealign.output.translations,
        reference = config["files"]["alignment_reference"],
        genemap = config["files"]["annotation"]
    output:
        mutation_summary = "pre-processed/{origin}/mutation_summary.tsv"
    log:
        "logs/mutation_summary_{origin}.txt"
    benchmark:
        "benchmarks/mutation_summary_{origin}.txt"
    params:
        outdir = "pre-processed/{origin}/translations",
        basename = "seqs",
        genes=config["genes"],
    conda: config["conda_environment"]
    shell:
        """
	conda activate nextstrain;
        python3 scripts/mutation_summary.py \
            --alignment {input.alignment} \
            --insertions {input.insertions} \
            --directory {params.outdir} \
            --basename {params.basename} \
            --reference {input.reference} \
            --genes {params.genes:q} \
            --genemap {input.genemap} \
            --output {output.mutation_summary} 2>&1 | tee {log}
        """

