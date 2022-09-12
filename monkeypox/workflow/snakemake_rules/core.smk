"""
This part of the workflow expects input files

        sequences = "data/sequences.fasta",
        metadata =  "data/metadata.tsv",

and will produce output files as

        auspice_json = auspice_dir + "/monkeypox_{build_name}.json"

Parameter are expected to sit in the `config` data structure.
In addition, `build_dir` and `auspice_dir` need to be defined upstream.
"""


rule wrangle_metadata:
    input:
        metadata="data/metadata.tsv",
    output:
        metadata=build_dir + "/{build_name}/metadata.tsv",
    params:
        strain_id=lambda w: config.get("strain_id_field", "strain"),
    shell:
        """
        python3 scripts/wrangle_metadata.py --metadata {input.metadata} \
                    --strain-id {params.strain_id} \
                    --output {output.metadata}
        """


rule filter:
    input:
        sequences="data/sequences.fasta",
        metadata=build_dir + "/{build_name}/metadata.tsv",
        exclude=config["exclude"],
        specific_exclude="config/{build_name}/exclude_accessions.txt",
        include="config/{build_name}/include_accessions.txt",
    output:
        sequences=build_dir + "/{build_name}/filtered.fasta",
        log=build_dir + "/{build_name}/filtered.log",
    params:
        min_date= lambda w: config[w.build_name]["min_date"],
        min_length=config["min_length"],
        exclude_where=lambda w: config[w.build_name]["exclude_where"],
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns strain \
            --exclude {input.exclude} {input.specific_exclude} \
            {params.exclude_where} \
            --include {input.include} \
            --output {output.sequences} \
            --min-length {params.min_length} \
            --min-date {params.min_date} \
            --output-log {output.log}
        """


rule align:
    message:
        """
        Aligning sequences to {input.reference}
        - filling gaps with N
        """
    input:
        sequences=rules.filter.output.sequences,
        reference="config/{build_name}/reference.fasta",
    output:
        alignment=build_dir + "/{build_name}/aligned.fasta",
        insertions=build_dir + "/{build_name}/insertions.fasta",
    params:
        max_indel=config["max_indel"],
        # seed_spacing = config["seed_spacing"]
        seed_spacing=500,
        terminal_bandwidth=500,
        excess_bandwidth=20,
    threads: workflow.cores
    shell:
        """
        nextalign run \
            --jobs {threads} \
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
            --include-reference \
            --output-insertions {output.insertions}
        """


rule mask:
    input:
        sequences=build_dir + "/{build_name}/aligned.fasta",
        mask="config/{build_name}/mask.bed",
    output:
        build_dir + "/{build_name}/masked.fasta",
    shell:
        """
        augur mask \
          --sequences {input.sequences} \
          --mask {input.mask} \
          --output {output}
        """


rule tree:
    message:
        "Building tree"
    input:
        alignment=build_dir + "/{build_name}/masked.fasta",
    output:
        tree=build_dir + "/{build_name}/tree_raw.nwk",
    threads: 3
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads auto \
            --tree-builder-args '-czb -redo -nt AUTO'
        """
    
rule fix_tree:
    message:
        "Building tree"
    input:
        tree=rules.tree.output.tree,
        alignment=build_dir + "/{build_name}/masked.fasta",
    output:
        tree=build_dir + "/{build_name}/tree_fixed.nwk",
    shell:
        """
        python3 scripts/fix_tree.py \
            --alignment {input.alignment} \
            --input-tree {input.tree} \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree
        """
    input:
        tree=rules.fix_tree.output.tree,
        alignment=build_dir + "/{build_name}/masked.fasta",
        metadata=build_dir + "/{build_name}/metadata.tsv",
    output:
        tree=build_dir + "/{build_name}/tree.nwk",
        node_data=build_dir + "/{build_name}/branch_lengths.json",
    params:
        root= lambda w: config[w.build_name]["root"],
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --root {params.root} \
            --divergence-unit mutations \
            --keep-polytomies \
            --output-node-data {output.node_data}
        """


rule ancestral:
    message:
        "Reconstructing ancestral sequences and mutations"
    input:
        tree=rules.refine.output.tree,
        alignment=build_dir + "/{build_name}/aligned.fasta",
    output:
        node_data=build_dir + "/{build_name}/nt_muts.json",
    params:
        inference="joint",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """


rule translate:
    message:
        "Translating amino acid sequences"
    input:
        tree=rules.refine.output.tree,
        node_data=rules.ancestral.output.node_data,
        genemap=config["genemap"],
    output:
        node_data=build_dir + "/{build_name}/aa_muts.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.genemap} \
            --output {output.node_data}
        """


rule clades:
    message:
        "Adding internal clade labels"
    input:
        tree=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        nuc_muts=rules.ancestral.output.node_data,
        clades="config/{build_name}/clades.tsv",  #customize
    output:
        node_data=build_dir + "/{build_name}/clades_raw.json",
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """


rule rename_clades:
    input:
        rules.clades.output.node_data,
    output:
        node_data=build_dir + "/{build_name}/clades.json",
    shell:
        """
        python scripts/clades_renaming.py \
        --input-node-data {input} \
        --output-node-data {output.node_data}
        """


rule colors:
    input:
        ordering="config/color_ordering.tsv",
        color_schemes="config/color_schemes.tsv",
        metadata=build_dir + "/{build_name}/metadata.tsv",
    output:
        colors=build_dir + "/{build_name}/colors.tsv",
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata} 2>&1
        """


rule export:
    message:
        "Exporting data files for for auspice"
    input:
        colors=rules.colors.output.colors,
        tree=rules.refine.output.tree,
        metadata=build_dir + "/{build_name}/metadata.tsv",
        branch_lengths=build_dir + "/{build_name}/branch_lengths.json",
        clades=rules.rename_clades.output.node_data,
        nt_muts=rules.ancestral.output.node_data,
        aa_muts=rules.translate.output.node_data,
        # description=config["description"],  #customize
        auspice_config="config/{build_name}/auspice_config.json",  #customize
    output:
        auspice_json=build_dir + "/{build_name}/raw_tree.json",
        root_sequence=build_dir + "/{build_name}/raw_tree_root-sequence.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence \
            --output {output.auspice_json}
        """


rule final_strain_name:
    input:
        auspice_json=build_dir + "/{build_name}/raw_tree.json",
        metadata=build_dir + "/{build_name}/metadata.tsv",
        root_sequence=build_dir + "/{build_name}/raw_tree_root-sequence.json",
    output:
        auspice_json=auspice_dir + "/{build_name}.json",
        root_sequence=auspice_dir + "/{build_name}_root-sequence.json",
    params:
        display_strain_field=lambda w: config.get("display_strain_field", "strain"),
    shell:
        """
        python3 scripts/set_final_strain_name.py --metadata {input.metadata} \
                --input-auspice-json {input.auspice_json} \
                --display-strain-name {params.display_strain_field} \
                --output {output.auspice_json}
        cp {input.root_sequence} {output.root_sequence}
        """
