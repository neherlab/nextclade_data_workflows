"""
This part of the workflow starts from files

  - pre-processed/sequences.fasta
  - pre-processed/metadata.tsv

and produces files

  - builds/{build_name}/sequences.fasta
  - builds/{build_name}/metadata.tsv

"""

build_dir = "builds"


rule synthetic_pick:
    input:
        counts="defaults/nr.tsv",
        metadata="pre-processed/open_pango_metadata.tsv",
    output:
        strains=build_dir + "/{build_name}/chosen_synthetic_strains.txt",
    log:
        "logs/synthetic_pick_{build_name}.txt",
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
        sequences="pre-processed/synthetic.fasta",
        strains=rules.synthetic_pick.output.strains,
    output:
        sequences=build_dir + "/{build_name}/sequences_raw.fasta",
    log:
        "logs/synthetic_select_{build_name}.txt",
    shell:
        """
        seqkit grep -f {input.strains} -o {output.sequences} {input.sequences} \
        2>&1 | tee {log}
        """


rule extract_metadata:
    input:
        strains=[
            build_dir + "/{build_name}/chosen_synthetic_strains.txt",
        ],
        metadata="data/metadata.tsv",
    output:
        metadata=build_dir + "/{build_name}/extracted_metadata.tsv",
    params:
        adjust=lambda w: config["builds"][w.build_name].get("metadata_adjustments", {}),
    benchmark:
        "benchmarks/extract_metadata_{build_name}.txt"
    run:
        import pandas as pd

        strains = set()
        for f in input.strains:
            with open(f) as fh:
                strains.update([x.strip() for x in fh if x[0] != "#"])

        d = pd.read_csv(input.metadata, index_col="strain", sep="\t")
        d = d[d.index.isin(list(strains))]
        if len(params.adjust):
            for adjustment in params.adjust:
                ind = d.eval(adjustment["query"])
                d.loc[ind, adjustment["dst"]] = d.loc[ind, adjustment["src"]]

        d.to_csv(output.metadata, sep="\t")


rule add_synthetic_metadata:
    input:
        metadata=rules.extract_metadata.output.metadata,
        synthetic=rules.synthetic_pick.output.strains,
    output:
        metadata=build_dir + "/{build_name}/metadata.tsv",
    log:
        "logs/add_synthetic_metadata_{build_name}.txt",
    shell:
        """
        python3 scripts/add_synthetic_metadata.py \
            --metadata {input.metadata} \
            --synthetic {input.synthetic} \
            --outfile {output.metadata} 2>&1 \
        | tee {log}
        """


rule exclude_outliers:
    input:
        sequences="builds/{build_name}/sequences_raw.fasta",
        metadata=rules.extract_metadata.output.metadata,
        exclude="profiles/exclude.txt",
    output:
        sampled_sequences="builds/{build_name}/sequences.fasta",
        sampled_strains="builds/{build_name}/subsample.txt",
    log:
        "logs/exclude_outliers_{build_name}.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sampled_sequences} \
            --output-strains {output.sampled_strains} \
        2>&1 | tee {log}
        """
