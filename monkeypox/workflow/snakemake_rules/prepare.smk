rule download:
    message:
        "Downloading sequences and metadata from data.nextstrain.org"
    output:
        sequences="data/sequences.fasta.xz",
        metadata="data/metadata.tsv.gz",
    params:
        sequences_url="https://data.nextstrain.org/files/workflows/monkeypox/sequences.fasta.xz",
        metadata_url="https://data.nextstrain.org/files/workflows/monkeypox/metadata.tsv.gz",
    shell:
        """
        curl -fsSL --compressed {params.sequences_url:q} --output {output.sequences}
        curl -fsSL --compressed {params.metadata_url:q} --output {output.metadata}
        """


rule decompress:
    message:
        "Decompressing sequences and metadata"
    input:
        sequences="data/sequences.fasta.xz",
        metadata="data/metadata.tsv.gz",
    output:
        sequences="data/sequences.fasta",
        metadata="data/metadata.tsv",
    shell:
        """
        gzip --decompress --keep {input.metadata}
        xz --decompress --keep {input.sequences}
        """


rule prealign:
    message:
        """
        Aligning sequences to {input.reference}
        - filling gaps with N
        """
    input:
        sequences=rules.decompress.output.sequences,
        reference="config/b1/reference.fasta",
    output:
        alignment="data/prealigned.fasta",
    params:
        max_indel=config["max_indel"],
        # seed_spacing = config["seed_spacing"]
        seed_spacing=500,
        terminal_bandwidth=500,
        excess_bandwidth=20,
    shell:
        """
        nextalign run \
            {input.sequences} \
            --retry-reverse-complement \
            --reference {input.reference} \
            --max-indel {params.max_indel} \
            --seed-spacing {params.seed_spacing} \
            --terminal-bandwidth {params.terminal_bandwidth} \
            --excess-bandwidth {params.excess_bandwidth} \
            --gap-alignment-side left \
            --genemap config/genemap.gff \
            --output-fasta {output.alignment} \
            --include-reference
        """


rule premask:
    input:
        sequences="data/prealigned.fasta",
        mask="config/b1/mask.bed",
    output:
        "data/premasked.fasta",
    shell:
        """
        augur mask \
            --sequences {input.sequences} \
            --mask {input.mask} \
            --output {output}
        """


rule deduplicate:
    """
    Remove identical sequences (even if they have differing Ns)
    Keep those sequences with fewer Ns
    """
    input:
        sequences=rules.premask.output,
    output:
        "data/duplicates.txt",
    shell:
        """
        python3 scripts/deduplicate.py \
            {input.sequences} \
            {output}
        """
