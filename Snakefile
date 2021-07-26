# rule preprocess:
#     input:
#         sequences = "pre-processed/vic.aligned.fasta",
#         metadata = "pre-processed/metadata.tsv",
#         sequence_index = "pre-processed/sequence_index.tsv",
#         mutation_summary = "pre-processed/mutation_summary.tsv"

rule download_clades:
    message: "Downloading clade definitions from {params.source} -> {output}"
    output:
        "data/clades.tsv"
    params:
        source = config["urls"]["clades"]
    shell: "curl {params.source} -o {output}"

rule download:
    output: temp("data/download_{year}.fasta")
    log: 'logs/download_{year}'
    shell: 
        "curl -ksSL -o {output} "
        "'https://www.viprbrc.org/brc/api/sequence?"
        "datatype=genome&"
        "completeseq=Y&"
        "family=influenza&"
        "flutype=B&"
        "fromyear={wildcards.year}&"
        "toyear={wildcards.year}&"
        "segment=4&"
        "host=human&"
        "metadata=ncbiAcc,continent,country,date,fluSeason,fluType,host,length,state,strainName&"
        "output=fasta' > {log}"

rule clean_download:
    input: rules.download.output
    output: "data/cleaned_download_{year}.fasta"
    shell: 
        """
        if [ $(wc -l < {input}) -lt 2 ]
        then
            touch {output}
        else
            cp {input} {output}
        fi
        """

rule join_downloads:
    input: 
        expand("data/cleaned_download_{year}.fasta",year=range(1980,2022))

    output: "data/download.fasta"
    shell: "cat {input} >> {output}"

rule parse:
    input:
        sequences = rules.join_downloads.output
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
        sequences = rules.parse.output.sequences,
        genemap = lambda w: config["files"]["annotation"][w.reference],
        reference = lambda w: config["files"]["alignment_reference"][w.reference]
    output:
        alignment = "pre-processed/{reference}.aligned.fasta",
        insertions = "pre-processed/{reference}.insertions.csv",
        translations = "pre-processed/{reference}.gene.HA.fasta"
    params:
        outdir = "pre-processed",
        genes = "HA",
    log:
        "logs/{reference}_prealign.txt"
    benchmark:
        "benchmarks/{reference}_align.txt"
    shell:
        """
        nextalign \
            --reference {input.reference} \
            --genemap {input.genemap} \
            --genes {params.genes} \
            --sequences {input.sequences} \
            --output-dir {params.outdir} \
            --output-basename {wildcards.reference}
        """

rule mutation_summary:
    message: "Summarizing {input.alignment}"
    input:
        alignment = rules.prealign.output.alignment,
        insertions = rules.prealign.output.insertions,
        translations = rules.prealign.output.translations,
        reference = lambda w: config["files"]["alignment_reference"][w.reference],
        genemap = lambda w: config["files"]["annotation"][w.reference]
    output:
        mutation_summary = "pre-processed/{reference}.mutation_summary.tsv"
    log:
        "logs/{reference}_mutation_summary.txt"
    benchmark:
        "benchmarks/{reference}_mutation_summary.txt"
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

rule enrich_metadata:
    input:
        metadata = rules.parse.output.metadata,
        mutation_summary = expand(rules.mutation_summary.output.mutation_summary,reference=['yam','vic'])
    output:
        enriched_metadata = "pre-processed/metadata_enriched.tsv"
    log: "logs/metadata_enrichment.txt"
    shell:
        """
        python3 scripts/metadata_enrichment.py \
            2>&1 | tee {log}
        """

rule subsample:
    input:
        sequences = expand("pre-processed/{subtype}.aligned.fasta",subtype=config['subtype']),
        metadata = "pre-processed/metadata_enriched.tsv",
    output:
        sequences = expand("build/{subtype}_subsample.fasta",subtype=config['subtype']),
        strains = expand("build/{subtype}_subsample.txt",subtype=config['subtype']),
    log:
        "logs/subsample.txt"
    benchmark:
        "benchmarks/subsample.txt"
    params:
        filter_arguments = config["filter"],
        include = config['refine']['root'],
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include-where ncbiAcc={params.include} \
            {params.filter_arguments} \
            --output {output.sequences} \
            --output-strains {output.strains} \
            2>&1 | tee {log}
        """


rule tree:
    message: "Building tree"
    input:
        alignment = rules.subsample.output.sequences
    output:
        tree = "build/tree_raw.nwk"
    params:
        args = lambda w: config["tree"].get("tree-builder-args","") if "tree" in config else ""
    log:
        "logs/tree.txt"
    benchmark:
        "benchmarks/tree.txt"
    threads: 16 
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
        alignment = rules.subsample.output.sequences,
        metadata = "pre-processed/metadata_enriched.tsv"
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
        divergence_unit = config["refine"].get("divergence_unit", 'mutations'),
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
            --divergence-unit {params.divergence_unit} \
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
        alignment = rules.subsample.output.sequences,
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

rule aa_muts_explicit:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        translations = expand(rules.prealign.output.translations,reference="vic")
    output:
        node_data = "build/aa_muts_explicit.json",
    params:
        genes = config.get('genes', 'HA')
    log:
        "logs/aamuts.txt"
    benchmark:
        "benchmarks/aamuts.txt"
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024)
    shell:
        """
        python3 scripts/explicit_translation.py \
            --tree {input.tree} \
            --translations {input.translations:q} \
            --genes {params.genes} \
            --output {output.node_data} 2>&1 | tee {log}
        """

rule clades:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.aa_muts_explicit.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = rules.download_clades.output
    output:
        node_data = "build/clades.json"
    log:
        "logs/clades.txt"
    benchmark:
        "benchmarks/clades.txt"
    resources:
        # Memory use scales primarily with size of the node data.
        mem_mb=lambda wildcards, input: 3 * int(input.size / 1024 / 1024)
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

rule export:
    message: "Exporting data files for auspice"
    input:
        tree = rules.refine.output.tree,
        # tree = "build/pruned_tree.nwk",
        metadata = rules.enrich_metadata.output.enriched_metadata,
        node_data = 
            [rules.__dict__[rule].output.node_data for rule in ['ancestral','refine','aa_muts_explicit','clades']],
        auspice_config = config['files']['auspice_config']
    output:
        auspice_json = "auspice/auspice.json",
    log:
        "logs/export.txt"
    benchmark:
        "benchmarks/export.txt"
    params:
        fields = "continent country fluSeason strainName"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data}\
            --include-root-sequence \
            --auspice-config {input.auspice_config} \
            --color-by-metadata {params.fields} \
            --output {output.auspice_json} 2>&1 | tee {log};
        """