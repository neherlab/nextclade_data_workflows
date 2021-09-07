'''
This part of the workflow starts from files

  - pre-processed/sequences.fasta
  - pre-processed/metadata.tsv

and produces files

  - builds/{build_name}/sequences.fasta
  - builds/{build_name}/metadata.tsv

'''

build_dir = "builds"

rule prepare_build:
    input:
        sequences = build_dir + "/{build_name}/sequences.fasta",
        metadata = build_dir + "/{build_name}/metadata.tsv"

rule subsample:
    message:
        """
        Subsample all sequences by '{wildcards.subsample}' scheme for build '{wildcards.build_name}' with the following parameters:
        """
    input:
        sequences = "data/sequences.fasta.xz",
        metadata = "data/metadata.tsv",
        sequence_index = "pre-processed/sequence_index.tsv",
        problematic_exclude = "pre-processed/problematic_exclude.txt",
        include = config["files"]["include"],
    output:
        sequences = build_dir + "/{build_name}/sample-{subsample}.fasta",
        strains=build_dir + "/{build_name}/sample-{subsample}.txt",
    log:
        "logs/subsample_{build_name}_{subsample}.txt"
    benchmark:
        "benchmarks/subsample_{build_name}_{subsample}.txt"
    params:
        filter_arguments = lambda w: config["builds"][w.build_name]["subsamples"][w.subsample]['filters'],
        date = (datetime.date.today() + datetime.timedelta(days=1)).strftime("%Y-%m-%d"),
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --sequence-index {input.sequence_index} \
            --include {input.include} \
            --exclude {input.problematic_exclude} \
            {params.filter_arguments} \
            --max-date {params.date} \
            --output {output.sequences} \
            --output-strains {output.strains} 2>&1 | tee {log}
        """

rule combine_subsamples:
    # Similar to rule combine_input_metadata, this rule should only be run if multiple inputs are being used (i.e. multiple origins)
    message:
        """
        Combine and deduplicate aligned & filtered FASTAs from multiple origins in preparation for subsampling: {input}.
        """
    input:
        lambda w: [build_dir + f"/{w.build_name}/sample-{subsample}.fasta"
                   for subsample in config["builds"][w.build_name]["subsamples"]]
    output:
        build_dir + "/{build_name}/sequences_raw.fasta"
    benchmark:
        "benchmarks/combine_subsamples_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/combine-and-dedup-fastas.py --input {input} --output {output}
        """

rule exclude_outliers:
    input:
        sequences = "builds/{build_name}/sequences_raw.fasta",
        metadata = "data/metadata.tsv",
        exclude = "profiles/exclude.txt",
                sequence_index = "pre-processed/sequence_index.tsv",
    output:
        sampled_sequences = "builds/{build_name}/sequences.fasta",
        sampled_strains = "builds/{build_name}/subsample.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --sequence-index {input.sequence_index} \
            --exclude {input.exclude} \
            --output {output.sampled_sequences} \
            --output-strains {output.sampled_strains}
        """
