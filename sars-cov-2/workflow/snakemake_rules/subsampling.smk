rule designated_lineages:
    input:
        sequences="pre-processed/synthetic.fasta",
    output:
        lineages="builds/designated_lineages.txt",
    shell:
        """
        seqkit seq -n {input.sequences} >{output.lineages}
        """


rule synthetic_pick:
    input:
        counts="defaults/nr.tsv",
        lineages=rules.designated_lineages.output.lineages,
        alias_file="pre-processed/alias.json",
        excluded_recombinants="profiles/clades/{build_name}/excluded_recombinants.txt",
        excludes="profiles/clades/excludes.txt",
    output:
        strains="builds/{build_name}/chosen_synthetic_strains.txt",
    params:
        strain_set=lambda w: config["strainSet"][w.build_name],
    shell:
        """
        python scripts/synthetic_pick.py \
            --counts {input.counts} \
            --lineages {input.lineages} \
            --alias-file {input.alias_file} \
            --build-name {params.strain_set} \
            --excluded-recombinants {input.excluded_recombinants} \
            --excludes {input.excludes} \
            --outfile {output.strains}
        """


rule synthetic_select:
    input:
        sequences="pre-processed/synthetic.fasta",
        strains=rules.synthetic_pick.output.strains,
    output:
        sequences="builds/{build_name}/sequences.fasta",
    shell:
        """
        seqkit grep -f {input.strains} -o {output.sequences} {input.sequences}
        """


rule add_synthetic_metadata:
    input:
        synthetic=rules.synthetic_pick.output.strains,
    output:
        metadata="builds/{build_name}/metadata.tsv",
    shell:
        """
        python3 scripts/add_synthetic_metadata.py \
            --synthetic {input.synthetic} \
            --outfile {output.metadata}
        """


rule add_designation_recency:
    input:
        designation_dates=rules.download_designation_dates.output.designation_dates,
    output:
        metadata="builds/{build_name}/designation_dates.tsv",
    shell:
        """
        python3 scripts/add_designation_recency.py \
            {input.designation_dates} \
            {output}
        """


rule add_designation_date_to_meta:
    input:
        metadata="builds/{build_name}/metadata.tsv",
        designation_dates=rules.add_designation_recency.output.metadata,
    output:
        metadata="builds/{build_name}/metadata_with_designation_date.tsv",
    shell:
        """
        tsv-select -H -f 1 {input.metadata} | \
        tsv-join -H --filter-file {input.designation_dates} \
            --key-fields 1 \
            --append-fields 2,3 \
        > {output}
        """


rule get_strains:
    input:
        sequences="builds/{build_name}/sequences.fasta",
    output:
        strains="builds/{build_name}/strains.txt",
    shell:
        """
        seqkit seq -n {input.sequences} >{output.strains}
        """
