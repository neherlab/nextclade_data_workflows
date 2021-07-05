rule preprocess:
    input:
        sequences = "pre-processed/seqs.aligned.fasta",
        metadata = "pre-processed/metadata.tsv",
        sequence_index = "pre-processed/sequence_index.tsv",
        mutation_summary = "pre-processed/mutation_summary.tsv"

rule download:
    output: "data/download.fasta"
    shell: 
        "curl -kL -o {output} "
        "'https://www.viprbrc.org/brc/api/sequence?"
        "datatype=genome&"
        "completeseq=Y&"
        "family=influenza&"
        "flutype=B&"
        "fromyear=2010&"
        "continent=Asia&"
        "segment=4&"
        "host=human&"
        "metadata=ncbiAcc,continent,country,date,fluSeason,fluType,host,length,state,strainName&"
        "output=fasta'"

rule parse:
    input:
        sequences = rules.download.output
    output:
        metadata = "pre-processed/metadata.tsv",
        sequences = "data/sequences.fasta"
    params:
        fields = "ncbiAcc continent country date fluSeason fluType host length state strainName"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --fields {params.fields} \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences}
        """

rule prealign:
    input:
        sequences = "data/sequences.fasta",
        genemap = config["files"]["annotation"],
        reference = config["files"]["alignment_reference"]
    output:
        alignment = "pre-processed/seqs.aligned.fasta",
        insertions = "pre-processed/seqs.insertions.csv",
        translations = "pre-processed/seqs.gene.HA.fasta"
    params:
        outdir = "pre-processed",
        genes = "HA",
    log:
        "logs/prealign.txt"
    benchmark:
        "benchmarks/align.txt"
    shell:
        """
        nextalign \
            --reference {input.reference} \
            --genemap {input.genemap} \
            --genes {params.genes} \
            --sequences {input.sequences} \
            --output-dir {params.outdir} \
            --output-basename "seqs"
        """

rule index_sequences:
    message:
        """
        Index sequence composition for faster filtering.
        """
    input:
        sequences = rules.preprocess.input.sequences
    output:
        sequence_index = rules.preprocess.input.sequence_index
    log:
        "logs/index_sequences.txt"
    benchmark:
        "benchmarks/index_sequences.txt"
    shell:
        """
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
        mutation_summary = rules.preprocess.input.mutation_summary
    log:
        "logs/mutation_summary.txt"
    benchmark:
        "benchmarks/mutation_summary.txt"
    params:
        outdir = "pre-processed/translations",
        basename = "seqs",
        genes= "HA"
    shell:
        """
        python3 scripts/mutation_summary.py \
            --alignment {input.alignment} \
            --insertions {input.insertions} \
            --directory {params.outdir} \
            --basename {params.basename} \
            --reference {input.reference} \
            --genes {params.genes} \
            --genemap {input.genemap} \
            --output {output.mutation_summary} 2>&1 | tee {log}
        """

rule mask:
    message:
        """
        Mask bases in alignment {input.alignment}
          - masking {params.mask_arguments}
        """
    input:
        alignment = rules.prealign.output.alignment
    output:
        alignment = "build/masked.fasta"
    log:
        "logs/mask.txt"
    benchmark:
        "benchmarks/mask.txt"
    params:
        mask_arguments = config.get("mask","")
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            {params.mask_arguments} \
            --output {output.alignment} 2>&1 | tee {log}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.mask.output.alignment
    output:
        tree = "build/tree_raw.nwk"
    params:
        args = lambda w: config["tree"].get("tree-builder-args","") if "tree" in config else ""
    log:
        "logs/tree.txt"
    benchmark:
        "benchmarks/tree.txt"
    threads: 8
    resources:
        # Multiple sequence alignments can use up to 40 times their disk size in
        # memory, especially for larger alignments.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 40 * int(input.size / 1024 / 1024)
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args {params.args} \
            --output {output.tree} \
            --nthreads {threads} 2>&1 | tee {log}
        """

rule refine:
    message:
        """
        Refining tree
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.mask.output.alignment,
        metadata = "pre-processed/metadata.tsv"
    output:
        tree = "build/tree.nwk",
        node_data = "build/branch_lengths.json"
    log:
        "logs/refine.txt"
    benchmark:
        "benchmarks/refine.txt"
    threads: 1
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024)
    params:
        root = config["refine"]["root"],
        coalescent = config["refine"].get("coalescent", "opt"),
        divergence_unit = config["refine"].get("divergence_unit", 'mutations'),
        clock_filter_iqd = config["refine"].get("clock_filter_iqd", 4),
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
            --coalescent {params.coalescent} \
            --divergence-unit {params.divergence_unit} \
            --date-confidence \
            --no-covariance \
           2>&1 | tee {log}
        """

rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - inferring ambiguous mutations
        """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.mask.output.alignment
    output:
        node_data = "build/nt_muts.json"
    log:
        "logs/ancestral.txt"
    benchmark:
        "benchmarks/ancestrald.txt"
    params:
        inference = config["ancestral"]["inference"]
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024)
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --infer-ambiguous 2>&1 | tee {log}
        """

rule export:
    message: "Exporting data files for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = "pre-processed/metadata.tsv",
        node_data = [rules.ancestral.output.node_data,rules.refine.output.node_data]
    output:
        auspice_json = "auspice/auspice.json",
    log:
        "logs/export.txt"
    benchmark:
        "benchmarks/export.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --include-root-sequence \
            --output {output.auspice_json} 2>&1 | tee {log};
        """