"""
This part of the workflow downloads files from S3

  - "data/{origin}/sequences.fasta.xz"
  - "data/{origin}/metadata.tsv"

and produces

  - pre-processed/filtered.fasta.xz
  - pre-processed/metadata.tsv

"""

import os


localrules:
    join_meta_nextclade,
    download_clade_emergence_dates,
    download_pango_alias,
    download_sequences,
    download_metadata,
    download_exclude,
    download_clades,
    preprocess,
    download_color_ordering,
    download_curated_pango,


rule preprocess:
    input:
        sequences="data/sequences.fasta.zst",
        metadata="data/metadata.tsv",
        sequence_index="pre-processed/sequence_index.tsv",
        problematic_exclude="pre-processed/problematic_exclude.txt",
        synthetic="pre-processed/synthetic.fasta",
    params:
        slack_hook=config.get("slackHook", "google.com"),
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
    if input.endswith(".zst"):
        return "zstd -dcq"
    return "cat"


rule download_sequences:
    message:
        "Downloading sequences from {params.address} -> {output[0]}"
    params:
        address=lambda w: config["origins"]["sequences"],
    output:
        "data/sequences.fasta.xz",
    shell:
        "aws s3 cp {params.address} {output}"

rule change_compression_to_zstd:
    input:
        "data/sequences.fasta.xz",
    output:
        "data/sequences.fasta.zst",
    threads: 5,
    shell:
        "xzcat {input} | zstd -c -10 -T4  > {output}"

rule download_metadata:
    message:
        "Downloading metadata from {params.address} -> {output}"
    params:
        deflate=lambda w: _infer_decompression(config["origins"]["metadata"]),
        address=lambda w: config["origins"]["metadata"],
    output:
        metadata="data/metadata_raw.tsv",
    shell:
        "aws s3 cp {params.address} - | {params.deflate} {input} > {output:q}"


rule fix_metadata:
    input:
        "data/metadata_raw.tsv",
    output:
        "data/metadata_raw2.tsv",
    shell:
        """
        awk  -F'\t' 'BEGIN {{OFS = FS}} {{if ($NF=="?") $NF="-inf"; print}}' {input} >{output}
        """


rule download_exclude:
    message:
        "Downloading exclude from {params.source} -> {output}"
    output:
        "data/exclude.txt",
    params:
        source=lambda w: config["origins"]["exclude"],
    shell:
        "curl {params.source} -o {output}"


rule download_clades:
    message:
        "Downloading clade definitions from {params.source} -> {output}"
    output:
        "builds/clades.tsv",
    params:
        source=config["data_source"]["clades"],
    shell:
        "curl {params.source} -o {output}"


rule download_color_ordering:
    message:
        "Downloading clade definitions from {params.source} -> {output}"
    output:
        "builds/color_ordering.tsv",
    params:
        source=config["data_source"]["color_ordering"],
    shell:
        "curl {params.source} -o {output}"


rule download_lat_longs:
    message:
        "Downloading clade definitions from {params.source} -> {output}"
    output:
        "builds/lat_longs.tsv",
    params:
        source=config["data_source"]["lat_longs"],
    shell:
        "curl {params.source} -o {output}"


rule download_curated_pango:
    output:
        "pre-processed/pango_raw.csv",
    params:
        source=config["data_source"]["pango"],
    shell:
        "curl {params.source} -o {output}"


rule download_pango_alias:
    output:
        "pre-processed/alias.json",
    params:
        source=config["data_source"]["aliases"],
    shell:
        "curl {params.source} -o {output}"


rule download_clade_emergence_dates:
    output:
        "pre-processed/clade_emergence_dates.tsv",
    params:
        source=config["data_source"]["clade_emergence_dates"],
    shell:
        "curl {params.source} -o {output}"

rule diagnostic:
    message:
        "Scanning metadata {input.metadata} for problematic sequences. Removing sequences with >{params.clock_filter} deviation from the clock and with more than {params.snp_clusters}."
    input:
        metadata="data/metadata_raw2.tsv",
        clade_emergence_dates="pre-processed/clade_emergence_dates.tsv",
    output:
        to_exclude="pre-processed/problematic_exclude.txt",
        exclude_reasons="pre-processed/exclude_reasons.txt",
    params:
        clock_filter=12,
        clock_filter_recent=17,
        clock_filter_lower_limit=-10,
        snp_clusters=1,
        rare_mutations=45,
        clock_plus_rare=50,
    log:
        "logs/diagnostics.txt",
    benchmark:
        "benchmarks/diagnostics.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000,
    shell:
        """
        python3 scripts/diagnostic.py \
            --metadata {input.metadata} \
            --clade_emergence_dates {input.clade_emergence_dates} \
            --clock-filter {params.clock_filter} \
            --rare-mutations {params.rare_mutations} \
            --clock-plus-rare {params.clock_plus_rare} \
            --snp-clusters {params.snp_clusters} \
            --output-exclusion-list {output.to_exclude} \
            --output-exclusion-reasons {output.exclude_reasons} \
            2>&1 | tee {log}
        """


rule index_sequences:
    message:
        """
        Index sequence composition for faster filtering.
        """
    input:
        sequences="data/sequences.fasta.xz",
    output:
        sequence_index="pre-processed/sequence_index.tsv",
    log:
        "logs/index_sequences.txt",
    benchmark:
        "benchmarks/index_sequences.txt"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index} 2>&1 | tee {log}
        """


rule nextclade_strainnames:
    message:
        "Extract strain names using tsv-select"
    input:
        "data/metadata_raw2.tsv",
    output:
        "pre-processed/metadata_strainnames.tsv",
    shell:
        """
        tsv-select -H -f strain {input} >{output}
        """


rule pango_strain_rename:
    message:
        "Convert pango strain names to nextclade strain names"
    input:
        metadata_strainnames="pre-processed/metadata_strainnames.tsv",
        pango="pre-processed/pango_raw.csv",
    output:
        pango_designations="pre-processed/pango_designations_nextstrain_names.csv",
        pango_designated_strains="pre-processed/pango_designated_strains_nextstrain_names.txt",
    shell:
        """
        python3 scripts/pango_strain_rename.py \
        --metadata-strainnames {input.metadata_strainnames} \
        --pango-in {input.pango} \
        --pango-designations {output.pango_designations} \
        --pango-designated-strains {output.pango_designated_strains} \
        2>&1
        """


rule fix_pango_lineages:
    message:
        "Add new column to open_pango_metadata_raw.tsv by joining pango_raw.csv on field strain name"
    input:
        metadata="data/metadata_raw2.tsv",
        pango_designations="pre-processed/pango_designations_nextstrain_names.csv",
    output:
        metadata="data/metadata.tsv",
    shell:
        """
        python3 scripts/fix_open_pango_lineages.py \
        --metadata {input.metadata} \
        --designations {input.pango_designations} \
        --output {output.metadata} \
        2>&1
        """


# rule open_pango:
#     # include only sequences that are in pango.csv using augur filter
#     # could be sped up using seqkit grep
#     # only problem: metadata, but can probably do with tsv-join
#     # and do inner join, not selecting sequences/meta without matches on both sides
#     input:
#         sequences = "data/sequences.fasta.xz",
#         pango = "pre-processed/pango_designated_strains_nextstrain_names.txt",
#         sequence_index = "pre-processed/sequence_index.tsv",
#         metadata = "data/metadata.tsv",
#     output:
#         sequences = "pre-processed/open_pango.fasta.xz",
#         metadata = "pre-processed/open_pango_metadata.tsv",
#         strains = "pre-processed/open_pango_strains.txt",
#     params:
#         date = (datetime.date.today() + datetime.timedelta(days=1)).strftime("%Y-%m-%d"),
#     log:
#         "logs/open_pango.txt"
#     benchmark:
#         "benchmarks/open_pango.txt"
#     conda: config["conda_environment"]
#     shell:
#         """
#         augur filter \
#             --sequences {input.sequences} \
#             --sequence-index {input.sequence_index} \
#             --metadata {input.metadata} \
#             --exclude-all \
#             --max-date {params.date} \
#             --include {input.pango} \
#             --output-metadata {output.metadata} \
#             --output-strains {output.strains} \
#             --output-sequences {output.sequences} \
#             2>&1 | tee {log}
#         """


rule get_designated_sequences:
    input:
        sequences="data/sequences.fasta.zst",
        pango="pre-processed/pango_designated_strains_nextstrain_names.txt",
    output:
        sequences="pre-processed/open_pango.fasta.zst",
    log:
        "logs/get_designated_sequences.txt",
    benchmark:
        "benchmarks/get_designated_sequences.txt"
    threads: 7
    shell:
        """
        zstdcat -T2 {input.sequences} | \
        seqkit grep -f {input.pango} 2>{log} | \
        zstd -c -10 -T4  >{output.sequences}
        """


rule get_designated_strains:
    input:
        sequences=rules.get_designated_sequences.output.sequences,
    output:
        strains="pre-processed/open_pango_strains.txt",
    log:
        "logs/get_designated_strains.txt",
    benchmark:
        "benchmarks/get_designated_strains.txt"
    threads: 3
    shell:
        """
        zstdcat -T2 {input.sequences} | \
        seqkit seq -in -o {output.strains} 2>{log}
        """


rule get_designated_metadata:
    input:
        strains="pre-processed/open_pango_strains.txt",
        metadata="data/metadata.tsv",
    output:
        metadata="pre-processed/open_pango_metadata.tsv",
    log:
        "logs/get_designated_metadata.txt",
    benchmark:
        "benchmarks/get_designated_metadata.txt"
    shell:
        """
        tsv-join -H --filter-file {input.strains} --key-fields 1 {input.metadata} 2>&1 >{output.metadata} | tee {log}
        """


rule strains:
    input:
        "data/metadata.tsv",
    output:
        "pre-processed/strains.txt",
    shell:
        "awk -F'\t' '{{print $1}}' {input} > {output}"


rule priorities:
    input:
        strains="pre-processed/strains.txt",
    output:
        priorities="pre-processed/priority.tsv",
    log:
        "logs/priorities.txt",
    benchmark:
        "benchmarks/priorities.txt"
    shell:
        """
        python3 scripts/priority_hash.py \
            --input {input.strains} \
            --output {output.priorities} \
            --seed 0 \
            2>&1 | tee {log} 
        """


rule join_meta_nextclade:
    input:
        open_pango_metadata=rules.get_designated_metadata.output.metadata,
    output:
        "pre-processed/full_sequence_details.tsv",
    shell:
        """
        tsv-select -H -f strain,pango_designated,deletions,insertions,substitutions {input.open_pango_metadata} > meta_muts.tsv
        aws s3 cp s3://nextstrain-ncov-private/nextclade.tsv.gz - \
            | zcat -d \
            | tsv-select -H -f seqName,missing,alignmentStart,alignmentEnd \
            | tsv-join -H -f meta_muts.tsv -k 1 -a 2-5 \
            > {output}
        rm meta_muts.tsv
        """


rule lineage_stats:
    input:
        reference="references/MN908947/reference.fasta",
        meta=rules.join_meta_nextclade.output,
    output:
        outfile="pre-processed/pango_matrix.npz",
    log:
        "logs/lineage_stats.txt",
    shell:
        """
        python3 scripts/lineage_matrix.py \
            --ref {input.reference} \
            --meta {input.meta} \
            --out {output.outfile} \
            2>&1 | tee {log} 
        """


rule make_synthetic_pangos:
    input:
        reference="references/MN908947/reference.fasta",
        matrix=rules.lineage_stats.output.outfile,
        alias=rules.download_pango_alias.output,
        overwrites="profiles/clades/lineage_overwrite.tsv",
    output:
        outfile="pre-processed/synthetic.fasta",
    log:
        "logs/make_synthetic_pangos.txt",
    shell:
        """
        python3 scripts/create_synthetic.py \
            --ref {input.reference} \
            --matrix {input.matrix} \
            --alias {input.alias} \
            --overwrites {input.overwrites} \
            --out {output.outfile} \
            2>&1 | tee {log} 
        """


# Need to create fake metadata for synthetic pangos
# Maybe in the above script?
