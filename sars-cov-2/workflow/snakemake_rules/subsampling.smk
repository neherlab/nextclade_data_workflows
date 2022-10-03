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
        metadata="pre-processed/open_pango_metadata.tsv.zst",
    output:
        strains=build_dir + "/{build_name}/chosen_synthetic_strains.txt",
    log:
        "logs/synthetic_pick_{build_name}.txt",
    shell:
        """
        zstdcat {input.metadata} | \
        python scripts/pick_synthetic.py \
            --designations /dev/stdin \
            --counts {input.counts} \
            --outfile {output.strains} 2>&1 \
        | tee {log}
        """


rule synthetic_select:
    input:
        sequences="pre-processed/synthetic.fasta",
        strains=rules.synthetic_pick.output.strains,
    output:
        sequences=build_dir + "/{build_name}/sequences.fasta",
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
        metadata=build_dir + "/{build_name}/metadata.tsv",
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
        sequences=build_dir + "/{build_name}/sequences.fasta",
    output:
        strains=build_dir + "/{build_name}/strains.txt",
    shell:
        """
        seqkit seq -n {input.sequences} >{output.strains}
        """
