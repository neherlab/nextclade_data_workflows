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
        "fromyear=2018&"
        "continent=Europe&"
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