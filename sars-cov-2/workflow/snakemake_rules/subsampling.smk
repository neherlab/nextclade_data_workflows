"""
This part of the workflow starts from files

  - pre-processed/sequences.fasta
  - pre-processed/metadata.tsv

and produces files

  - builds/{build_name}/sequences.fasta
  - builds/{build_name}/metadata.tsv

"""


rule designated_lineages:
    input:
        sequences="pre-processed/synthetic.fasta",
    output:
        lineages="builds" + "/{build_name}/designated_lineages.txt",
    shell:
        """
        seqkit seq -n {input.sequences} >{output.lineages}
        """


rule synthetic_pick:
    input:
        counts="defaults/nr.tsv",
        lineages=rules.designated_lineages.output.lineages,
        alias_file="pre-processed/alias.json",
    output:
        strains="builds" + "/{build_name}/chosen_synthetic_strains.txt",
    shell:
        """
        python scripts/synthetic_pick.py \
            --counts {input.counts} \
            --lineages {input.lineages} \
            --alias-file {input.alias_file} \
            --outfile {output.strains}
        """


rule synthetic_select:
    input:
        sequences="pre-processed/synthetic.fasta",
        strains=rules.synthetic_pick.output.strains,
    output:
        sequences="builds" + "/{build_name}/sequences.fasta",
    log:
        "logs/synthetic_select_{build_name}.txt",
    shell:
        """
        seqkit grep -f {input.strains} -o {output.sequences} {input.sequences} \
        2>&1 | tee {log}
        """


rule add_synthetic_metadata:
    input:
        synthetic=rules.synthetic_pick.output.strains,
    output:
        metadata="builds" + "/{build_name}/metadata.tsv",
    log:
        "logs/add_synthetic_metadata_{build_name}.txt",
    shell:
        """
        python3 scripts/add_synthetic_metadata.py \
            --synthetic {input.synthetic} \
            --outfile {output.metadata} 2>&1 \
        | tee {log}
        """


rule get_strains:
    input:
        sequences="builds" + "/{build_name}/sequences.fasta",
    output:
        strains="builds" + "/{build_name}/strains.txt",
    shell:
        """
        seqkit seq -n {input.sequences} >{output.strains}
        """
