'''
This part of the workflow downloads reference data sets

  - "reference-datasets/{build_name}/sequences.fasta"
  - "reference-datasets/{build_name}/metadata.tsv"

and combines them with custom data to produce

  - builds/{build_name}/sequences.fasta
  - builds/{build_name}/metadata.tsv

'''
rule prepare_reference_build:
    input:
        sequences = "builds/{build_name}/sequences.fasta",
        metadata = "builds/{build_name}/metadata.tsv"


def _infer_decompression(input):
    """
    Returns a shell command to decompress the piped stream,
    which will itself produce a stream of decompressed data to stdout.
    If no decompression is needed, returns `cat`.
    NOTE: a lot of this will become unnecessary once `augur` handles
    compressed sequence inputs.
    """
    if input.endswith(".xz"):
        return "xz -dcq"
    if input.endswith(".gz"):
        return "gunzip -cq"
    return "cat"


rule download:
    message: "Downloading reference sequences and metadata"
    output:
        sequences = "reference-datasets/{build_name}/sequences.fasta",
        metadata =  "reference-datasets/{build_name}/metadata.fasta"
    params:
        meta = lambda w: config["builds"][w.build_name]["reference_metadata"],
        seq =  lambda w: config["builds"][w.build_name]["reference_sequences"],
        deflate_seq = lambda w: _infer_decompression(config["builds"][w.build_name]["reference_metadata"]),
        deflate_meta = lambda w: _infer_decompression(config["builds"][w.build_name]["reference_sequences"])
    conda: config["conda_environment"]
    shell:
        """
        curl {params.meta} | {params.deflate_meta} > {output.metadata:q}
        curl {params.seq} | {params.deflate_seq} > {output.sequences:q}
        """


rule combine_input_metadata:
    # this rule is intended to be run _only_ if we have defined multiple inputs ("origins")
    message:
        """
        Combining metadata files {input.ref_metadata} and {input.user_metadata} and adding columns to represent origin
        """
    input:
        ref_metadata = rules.download.output.metadata,
        user_metadata = lambda w: config["builds"][w.build_name]["user_metadata"]
    output:
        metadata = rules.prepare_reference_build.input.metadata
    log:
        "logs/combine_input_metadata_{build_name}.txt"
    benchmark:
        "benchmarks/combine_input_metadata_{build_name}.txt"
    run:
        import pandas as pd

        ref_meta = pd.read_csv(input.ref_metadata, sep='\t')
        ref_meta["source"] = "background"

        user_metadata_files = [input.user_metadata] if type(input.user_metadata)==str else input.user_metadata
        user_metadata = []
        for fname in user_metadata_files:
            tmp = pd.read_csv(fname, sep=None)
            tmp['source'] = "foreground"
            user_metadata.append(tmp)

        combined_meta = pd.concat([ref_meta]+user_metadata, axis=0).fillna("?").drop_duplicates(subset="strain", keep='first')
        combined_meta.to_csv(output.metadata, sep='\t')


rule combine_sequences:
    # Similar to rule combine_input_metadata, this rule should only be run if multiple inputs are being used (i.e. multiple origins)
    message:
        """
        Combine and deduplicate aligned & filtered FASTAs from multiple origins in preparation for subsampling: {input}.
        """
    input:
        lambda w: [f"reference-datasets/{w.build_name}/sequences.fasta"] + ([config["builds"][w.build_name]["user_sequences"]] if type(config["builds"][w.build_name]["user_sequences"])==str else config["builds"][w.build_name]["user_sequences"])
    output:
        rules.prepare_reference_build.input.sequences
    benchmark:
        "benchmarks/combine_sequences_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/combine-and-dedup-fastas.py --input {input} --output {output}
        """

