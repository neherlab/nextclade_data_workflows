# rule preprocess:
#     input:
#         sequences = "pre-processed/vic.aligned.fasta",
#         metadata = "pre-processed/metadata.tsv",
#         sequence_index = "pre-processed/sequence_index.tsv",
#         mutation_summary = "pre-processed/mutation_summary.tsv"
wildcard_constraints:
    flu_type="[AB]",
    segment="\d",
    year="\d\d\d\d",


rule download_clades:
    message: "Downloading clade definitions for {wildcards.strain} from {params.source} -> {output}"
    output:
        "data/clades_{strain}.tsv"
    params:
        source = lambda w: config["urls"]["clades"][w.strain]
    shell: "curl {params.source} -o {output}"

rule download:
    output: temp("data/download_{flu_type}_{segment}_{year}.fasta")
    log: 'logs/download_{flu_type}_{segment}_{year}'
    shell:
        "curl -ksSL -o {output} "
        "'https://www.viprbrc.org/brc/api/sequence?"
        "datatype=genome&"
        "completeseq=Y&"
        "family=influenza&"
        "flutype={wildcards.flu_type}&"
        "fromyear={wildcards.year}&"
        "toyear={wildcards.year}&"
        "segment={wildcards.segment}&"
        "host=human&"
        "metadata=ncbiAcc,continent,country,date,fluSeason,fluType,host,length,state,strainName&"
        "output=fasta' > {log}"

rule clean_download:
    input: rules.download.output
    output: temp("data/cleaned_download_{flu_type}_{segment}_{year}.fasta")
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
        expand("data/cleaned_download_{{flu_type}}_{{segment}}_{year}.fasta",year=range(1989,2022))
    output: "data/download_{flu_type}_{segment}.fasta"
    shell: "cat {input} >> {output}"

rule parse:
    input:
        sequences = rules.join_downloads.output
    output:
        metadata = "pre-processed/metadata_{flu_type}_{segment}.tsv",
        sequences = "data/sequences_{flu_type}_{segment}.fasta"
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

genes = ["SigPep","HA1","HA2"]

def flu_type(reference):
    if reference in ["vic","yam"]:
        return "B"
    return "A"

def genes(reference,segment):
    segment = int(segment)
    if reference in ["vic","yam"] and segment <= 2:
        if segment == 1:
            segment = 2
        else:
            segment = 1 
    genes = [
        ["PB2"],
        ["PB1"],
        ["PA"],
        ["SigPep","HA1","HA2"],
        ["NP"],
        ["NA"],
        ["MA"],
        ["MS"],
    ]
    seg_name = [
        "PB2",
        "PB1",
        "PA",
        "HA",
        "NP",
        "NA",
        "MA",
        "MS",
    ]
    return {'genes': genes[segment-1], 'seg_name': seg_name[segment-1].lower()}

rule prealign:
    input:
        sequences = lambda w: expand(rules.parse.output.sequences,flu_type=flu_type(w.reference),segment=w.segment),
        genemap = lambda w: f"references/{w.reference}/{genes(w.reference,w.segment)['seg_name']}/genemap.gff",
        reference = lambda w: f"references/{w.reference}/{genes(w.reference,w.segment)['seg_name']}/reference.fasta"
    output:
        alignment = "pre-processed/{reference}_{segment}.aligned.fasta",
        insertions = "pre-processed/{reference}_{segment}.insertions.csv",
        # translations = "pre-processed/{reference}_{segment}.gene.{params.genes}.fasta"
    params:
        outdir = "pre-processed",
        basename = lambda w: f"{w.reference}_{w.segment}",
        genes = lambda w: ",".join(genes(w.reference,w.segment)['genes']),
    log:
        "logs/{reference}_{segment}_prealign.txt"
    benchmark:
        "benchmarks/{reference}_{segment}_align.txt"
    shell:
        """
        nextalign \
            --reference {input.reference} \
            --genemap {input.genemap} \
            --genes {params.genes} \
            --sequences {input.sequences} \
            --output-dir {params.outdir} \
            --output-basename {params.basename} \
        > {log} 2>&1
        """

rule mutation_summary:
    message: "Summarizing {input.alignment}"
    input:
        alignment = rules.prealign.output.alignment,
        insertions = rules.prealign.output.insertions,
        # translations = lambda w: expand("pre-processed/{{reference}}_{{segment}}.gene.{genes}.fasta",genes=genes(w.segment)),
        genemap = lambda w: f"references/{w.reference}/{genes(w.reference,w.segment)['seg_name']}/genemap.gff",
        reference = lambda w: f"references/{w.reference}/{genes(w.reference,w.segment)['seg_name']}/reference.fasta"
    output:
        mutation_summary = "pre-processed/{reference}_{segment}.mutation_summary.tsv"
    log:
        "logs/{reference}_{segment}_mutation_summary.txt"
    benchmark:
        "benchmarks/{reference}_{segment}_mutation_summary.txt"
    params:
        outdir = "pre-processed",
        basename = lambda w: f"{w.reference}_{w.segment}",
        genes = lambda w: ",".join(genes(w.reference,w.segment)['genes']),
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
        metadata = lambda w:  expand(rules.parse.output.metadata,segment="4",flu_type=w.flu_type),
        mutation_summary = lambda w: expand(rules.mutation_summary.output.mutation_summary,reference=(['yam','vic'] if w.flu_type == 'B' else ['h1','h3']),segment="4",flu_type=w.flu_type)
    output:
        enriched_metadata = "pre-processed/metadata_enriched_{flu_type}.tsv"
    log: "logs/metadata_enrichment_{flu_type}.txt"
    shell:
        """
        python3 scripts/metadata_enrichment.py \
            --flu-type {wildcards.flu_type} \
            2>&1 | tee {log}
        """

rule subsample:
    input:
        sequences = "pre-processed/{strain}_{segment}.aligned.fasta",
        metadata = lambda w: f"pre-processed/metadata_enriched_{flu_type(w.strain)}.tsv",
    output:
        sequences = "build/{strain}/{segment}/subsample.fasta",
        strains = "build/{strain}/{segment}/subsample.txt",
    log:
        "logs/subsample_{strain}_{segment}.txt",
    params:
        filter_arguments = lambda w: config["filter"][w.strain],
        include = lambda w: config['refine']['root'][w.strain],
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
        tree = "build/{strain}/{segment}/tree_raw.nwk"
    params:
        args = lambda w: config["tree"].get("tree-builder-args","") if "tree" in config else ""
    log:
        "logs/tree_{strain}_{segment}.txt"
    threads: 1
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
        metadata = lambda w: f"pre-processed/metadata_enriched_{flu_type(w.strain)}.tsv"
    output:
        tree = "build/{strain}/{segment}/tree.nwk",
        node_data = "build/{strain}/{segment}/branch_lengths.json"
    log:
        "logs/refine_{strain}_{segment}.txt"
    threads: 1
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024)
    params:
        root = lambda w: config["refine"]["root"][w.strain],
        divergence_unit = "mutations",
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
        node_data = "build/{strain}/{segment}/nt_muts.json"
    log:
        "logs/ancestral_{strain}_{segment}.txt"
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
        # translations = lambda w: expand("pre-processed/{{reference}}_{{segment}}.gene.{genes}.fasta",genes=genes(w.segment)),
    output:
        node_data = "build/{strain}/{segment}/aa_muts_explicit.json",
        # translations = vic.gene.HA_withInternalNodes.fasta
    params:
        genes = lambda w: genes(w.strain,w.segment)['genes'],
        translations = lambda w: expand("pre-processed/{strain}_{segment}.gene.{genes}.fasta",strain=w.strain,segment=w.segment,genes=genes(w.strain,w.segment)['genes']),
    log:
        "logs/aamuts_{strain}_{segment}.txt"
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024)
    shell:
        """
        python3 scripts/explicit_translation.py \
            --tree {input.tree} \
            --translations {params.translations:q} \
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
        node_data = "build/{strain}/{segment}/clades.json"
    log:
        "logs/clades_{strain}_{segment}.txt"
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
        metadata = lambda w: expand(rules.enrich_metadata.output.enriched_metadata,flu_type=flu_type(w.strain)),
        node_data = 
            [rules.__dict__[rule].output.node_data for rule in ['ancestral','refine','aa_muts_explicit','clades']],
        auspice_config = config['files']['auspice_config']
    output:
        auspice_json = "auspice/{strain}/{segment}/auspice.json",
    log:
        "logs/export_{strain}_{segment}.txt"
    params:
        fields = "continent country fluSeason strainName"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=100
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