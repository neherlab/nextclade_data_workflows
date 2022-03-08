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
        priority = "pre-processed/priority.tsv",
    output:
        sequences = build_dir + "/{build_name}/sample-{subsample,[^-]*}.fasta",
        strains=build_dir + "/{build_name}/sample-{subsample,[^-]*}.txt",
    log:
        "logs/subsample_{build_name}_{subsample}.txt"
    benchmark:
        "benchmarks/subsample_{build_name}_{subsample}.txt"
    params:
        filter_arguments = lambda w: config["builds"][w.build_name]["subsamples"][w.subsample]['filters'],
        date = (datetime.date.today() + datetime.timedelta(days=1)).strftime("%Y-%m-%d"),
        exclude_where_args = config["exclude-where-args"],
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
            {params.filter_arguments} {params.exclude_where_args} \
            --query "pango_lineage != ''" \
            --priority {input.priority} \
            --output {output.sequences} \
            --output-strains {output.strains} 2>&1 | tee {log}
        """

# 1. Choose synthetic sequences
# 2. Choose non-synthetic sequences


rule pango_pick:
    input:
        counts = "defaults/nr.tsv",
        metadata = "pre-processed/open_pango_metadata.tsv",
        exclude = "pre-processed/problematic_exclude.txt",
    output:
        strains = build_dir + "/{build_name}/chosen_pango_strains.txt",
    log:
        "logs/pango_pick_{build_name}.txt"
    shell:
        """
        python scripts/pick_samples.py \
            --designations {input.metadata} \
            --counts {input.counts} \
            --exclude {input.exclude} \
            --outfile {output.strains} 2>&1 \
        | tee {log}
        """

rule pango_select:
    input:
        sequences = "pre-processed/open_pango.fasta.xz",
        strains = rules.pango_pick.output.strains,
    output:
        sequences = build_dir + "/{build_name}/picked_pango.fasta",
    shell:
        """
        xzcat {input.sequences} | \
        seqkit grep -f {input.strains} -o {output.sequences}
        """

rule pango_sampling:
    input:
        sequences = rules.pango_select.output.sequences,
        sequence_index = "pre-processed/sequence_index.tsv",
        metadata = "pre-processed/open_pango_metadata.tsv",
    output:
        sequences = build_dir + "/{build_name}/sample-{subsample}-pango.fasta",
        strains=build_dir + "/{build_name}/sample-{subsample}-pango.txt",
    log:
        "logs/subsample_{build_name}_{subsample}-pango.txt"
    benchmark:
        "benchmarks/subsample_{build_name}_{subsample}-pango.txt"
    params:
        exclude_where_args = config["exclude-where-args"],
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --exclude-where Nextstrain_clade='21K (Omicron)' Nextstrain_clade='21L (Omicron)' Nextstrain_clade='21M (Omicron)' recombinant=True \
            --output {output.sequences} \
            --output-strains {output.strains} 2>&1 | tee {log}
        """
    

rule synthetic_pick:
    input:
        counts = "defaults/nr.tsv",
        metadata = "pre-processed/open_pango_metadata.tsv",
    output:
        strains = build_dir + "/{build_name}/chosen_synthetic_strains.txt",
    log:
        "logs/synthetic_pick_{build_name}.txt"
    shell:
        """
        python scripts/pick_synthetic.py \
            --designations {input.metadata} \
            --counts {input.counts} \
            --outfile {output.strains} 2>&1 \
        | tee {log}
        """

rule synthetic_select:
    input:
        sequences = "pre-processed/synthetic.fasta",
        strains = rules.synthetic_pick.output.strains,
    output:
        sequences = build_dir + "/{build_name}/picked_synthetic.fasta",
    log:
        "logs/synthetic_select_{build_name}.txt"
    shell:
        """
        seqkit grep -f {input.strains} -o {output.sequences} {input.sequences} \
        2>&1 | tee {log}
        """


rule combine_subsamples:
    # Similar to rule combine_input_metadata, this rule should only be run if multiple inputs are being used (i.e. multiple origins)
    message:
        """
        Combine and deduplicate aligned & filtered FASTAs from multiple origins in preparation for subsampling: {input}.
        """
    input:
        natural = lambda w: [build_dir + f"/{w.build_name}/sample-{subsample}.fasta"
                   for subsample in config["builds"][w.build_name]["subsamples"]],
        synthetic = rules.synthetic_select.output.sequences,
    output:
        build_dir + "/{build_name}/sequences_raw.fasta"
    benchmark:
        "benchmarks/combine_subsamples_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/combine-and-dedup-fastas.py --input {input} --output {output}
        """

rule extract_metadata:
    input:
        strains = lambda w: [build_dir + f"/{w.build_name}/sample-{subsample}.txt"
                for subsample in config["builds"][w.build_name]["subsamples"]],
        metadata = "data/metadata.tsv"
    output:
        metadata = build_dir + "/{build_name}/extracted_metadata.tsv"
    params:
        adjust = lambda w: config["builds"][w.build_name].get("metadata_adjustments",{}),
    benchmark:
        "benchmarks/extract_metadata_{build_name}.txt"
    run:
        import pandas as pd
        strains = set()
        for f in input.strains:
            with open(f) as fh:
                strains.update([x.strip() for x in fh if x[0]!='#'])

        d = pd.read_csv(input.metadata, index_col='strain', sep='\t')
        d = d[d.index.isin(list(strains))]
        if len(params.adjust):
            for adjustment  in params.adjust:
                ind = d.eval(adjustment["query"])
                d.loc[ind, adjustment['dst']] = d.loc[ind, adjustment['src']]

        d.to_csv(output.metadata, sep='\t')

rule add_synthetic_metadata:
    input:
        metadata = rules.extract_metadata.output.metadata,
        synthetic = rules.synthetic_pick.output.strains,
    output:
        metadata = rules.prepare_build.input.metadata
    log:
        "logs/add_synthetic_metadata_{build_name}.txt"
    shell:
        """
        python3 scripts/add_synthetic_metadata.py \
            --metadata {input.metadata} \
            --sequences {input.synthetic} \
            --outfile {output.metadata} 2>&1 \
        | tee {log}
        """

rule exclude_outliers:
    input:
        sequences = "builds/{build_name}/sequences_raw.fasta",
        metadata = rules.prepare_build.input.metadata,
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