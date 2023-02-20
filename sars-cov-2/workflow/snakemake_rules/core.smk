localrules:
    add_branch_labels,
    colors,
    internal_pango,
    overwrite_recombinant_clades,
    add_recombinants_to_tree,
    remove_recombinants_from_alignment,
    identify_recombinants,


build_dir = config.get("build_dir", "builds")
auspice_dir = config.get("auspice_dir", "auspice")
auspice_prefix = config.get("auspice_prefix", "ncov")


rule align:
    input:
        sequences="builds/nextclade/sequences.fasta",
        genemap=config["files"]["annotation"],
        reference=config["files"]["alignment_reference"],
    output:
        alignment="builds/nextclade/aligned.fasta",
        translations=expand(
            "builds/nextclade/translations/aligned.gene.{gene}.fasta",
            gene=config.get("genes", ["S"]),
        ),
    params:
        outdir=lambda w: "builds/nextclade/translations/aligned.gene.{gene}.fasta",
        genes="ORF1a,ORF1b,S,ORF3a,M,N",
        basename="aligned",
    threads: 4
    shell:
        """
        nextalign run \
            --jobs={threads} \
            --input-ref {input.reference} \
            --input-gene-map {input.genemap} \
            --genes {params.genes} \
            {input.sequences} \
            --output-translations {params.outdir} \
            --output-fasta {output.alignment} \
            2>&1
        """


rule mask:
    input:
        alignment=rules.align.output.alignment,
    output:
        alignment="builds/nextclade/masked.fasta",
    params:
        mask_arguments=lambda w: config.get("mask", ""),
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            {params.mask_arguments} \
            --output {output.alignment} 2>&1
        """


rule separate_recombinants:
    input:
        alignment=rules.mask.output.alignment,
        alias_json=rules.download_pango_alias.output,
    output:
        without_recombinants="builds/nextclade/masked_without_recombinants.fasta",
        recombinant_alignments=expand(
            "builds/nextclade/masked_recombinant_{tree_recombinants}.fasta",
            tree_recombinants=config["tree-recombinants"],
        ),
        recombinants="builds/nextclade/recombinants.txt",
    params:
        tree_recombinants=",".join(config["tree-recombinants"]),
    shell:
        """
        python3 scripts/separate_recombinants.py \
            --alignment {input.alignment} \
            --output-without-recombinants {output.without_recombinants} \
            --alias-json {input.alias_json} \
            --tree-recombinants {params.tree_recombinants} \
            --recombinants {output.recombinants}
        """


rule tree:
    input:
        alignment="builds/nextclade/masked_without_recombinants.fasta",
        constraint_tree=config["files"]["constraint_tree"],
        exclude_sites=config["files"]["exclude_sites"],
    output:
        tree="builds/nextclade/tree_raw.nwk",
    params:
        args="'-czb -g defaults/constraint.nwk'",
    threads: 8
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --exclude-sites {input.exclude_sites} \
            --tree-builder-args {params.args} \
            --output {output.tree} \
            --nthreads {threads} 2>&1
        """


rule recombinant_tree:
    input:
        alignment="builds/nextclade/masked_recombinant_{recombinant}.fasta",
    output:
        tree="builds/nextclade/tree_raw_recombinant_{recombinant}.nwk",
    params:
        constraint=lambda w: f"-g profiles/clades/constraint_{w.recombinant}.nwk"
        if os.path.exists(f"profiles/clades/constraint_{w.recombinant}.nwk")
        else "",
    shell:
        """
        # Check if there are 3 or more sequences in the alignment
        if [ $(grep -c ">" {input.alignment}) -lt 3 ]; then
            python scripts/simple_tree.py \
                --alignment {input.alignment} \
                --tree {output.tree} 2>&1
        else 
            augur tree \
                --alignment {input.alignment} \
                --tree-builder-args "-czb {params.constraint}" \
                --output {output.tree} \
                2>&1
        fi
        """


rule add_recombinants_to_tree:
    """
    Adding recombinants to root of raw tree
    """
    input:
        tree=rules.tree.output.tree,
        recombinants="builds/nextclade/recombinants.txt",
        recombinant_trees=expand(
            "builds/nextclade/tree_raw_recombinant_{recombinant}.nwk",
            recombinant=config["tree-recombinants"],
        ),
    output:
        tree="builds/nextclade/tree_with_recombinants.nwk",
    params:
        root=config["refine"]["root"],
        joined_trees=lambda w: ",".join(
            rules.add_recombinants_to_tree.input.recombinant_trees
        ),
    shell:
        """
        python scripts/add_recombinants.py \
            --tree {input.tree} \
            --recombinants {input.recombinants} \
            --recombinant-trees {params.joined_trees} \
            --root {params.root} \
            --output {output.tree} 2>&1 | tee {log}
        """


rule refine:
    input:
        tree=rules.add_recombinants_to_tree.output.tree,
        alignment=rules.align.output.alignment,
        metadata="builds/{build_name}/metadata.tsv",
    output:
        tree=build_dir + "/{build_name}/tree.nwk",
        node_data=build_dir + "/{build_name}/branch_lengths.json",
    log:
        "logs/refine_{build_name}.txt",
    benchmark:
        "benchmarks/refine_{build_name}.txt"
    params:
        root=config["refine"]["root"],
        divergence_unit=config["refine"].get("divergence_unit", "mutations"),
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --keep-root \
            --divergence-unit {params.divergence_unit} | tee {log}
        """


rule ancestral:
    input:
        tree=rules.refine.output.tree,
        alignment=rules.align.output.alignment,
    output:
        node_data=build_dir + "/{build_name}/nt_muts.json",
    log:
        "logs/ancestral_{build_name}.txt",
    benchmark:
        "benchmarks/ancestral_{build_name}.txt"
    params:
        inference="joint",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --infer-ambiguous 2>&1 | tee {log}
        """


rule translate:
    input:
        tree=rules.refine.output.tree,
        node_data=rules.ancestral.output.node_data,
        reference=config["files"]["reference"],
    output:
        node_data=build_dir + "/{build_name}/aa_muts.json",
    log:
        "logs/translate_{build_name}.txt",
    benchmark:
        "benchmarks/translate_{build_name}.txt"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """


rule aa_muts_explicit:
    input:
        tree=rules.refine.output.tree,
        translations=lambda w: rules.align.output.translations,
    output:
        node_data=build_dir + "/{build_name}/aa_muts_explicit.json",
        translations=expand(
            build_dir
            + "/{{build_name}}/translations/aligned.gene.{gene}_withInternalNodes.fasta",
            gene=config.get("genes", ["S"]),
        ),
    params:
        genes=config.get("genes", "S"),
    log:
        "logs/aamuts_{build_name}.txt",
    benchmark:
        "benchmarks/aamuts_{build_name}.txt"
    shell:
        """
        python3 scripts/explicit_translation.py \
            --tree {input.tree} \
            --translations {input.translations:q} \
            --genes {params.genes} \
            --output {output.node_data} 2>&1 | tee {log}
        """


rule internal_pango:
    input:
        tree=rules.refine.output.tree,
        alias=rules.download_pango_alias.output,
        synthetic=rules.synthetic_pick.output,
        designations=rules.pango_strain_rename.output.pango_designations,
    output:
        node_data=build_dir + "/{build_name}/internal_pango.json",
    log:
        "logs/internal_pango_{build_name}.txt",
    shell:
        """
        python scripts/internal_pango.py \
            --tree {input.tree} \
            --synthetic {input.synthetic} \
            --alias {input.alias} \
            --designations {input.designations} \
            --output {output.node_data} \
            --field-name Nextclade_pango 2>&1 | tee {log}
        """


rule clades_legacy:
    input:
        tree=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        nuc_muts=rules.ancestral.output.node_data,
        clades="builds/clades.tsv",
        internal_pango=rules.internal_pango.output.node_data,
        alias=rules.download_pango_alias.output,
    output:
        node_data=build_dir + "/{build_name}/clades_legacy.json",
    log:
        "logs/clades_{build_name}.txt",
    benchmark:
        "benchmarks/clades_{build_name}.txt"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data clades_raw.tmp 2>&1 | tee {log}
        python scripts/overwrite_recombinant_clades.py \
            --clades clades_raw.tmp \
            --internal-pango {input.internal_pango} \
            --alias {input.alias} \
            --clade-type clade_legacy \
            --output {output.node_data}
        rm clades_raw.tmp
        sed -i'' 's/clade_membership/clade_legacy/gi' {output.node_data}
        """


rule clades:
    input:
        tree=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        nuc_muts=rules.ancestral.output.node_data,
        clades="builds/clades_nextstrain.tsv",
        internal_pango=rules.internal_pango.output.node_data,
        alias=rules.download_pango_alias.output,
    output:
        node_data=build_dir + "/{build_name}/clades.json",
        node_data_nextstrain=build_dir + "/{build_name}/clades_nextstrain.json",
    log:
        "logs/clades_{build_name}.txt",
    benchmark:
        "benchmarks/clades_{build_name}.txt"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data clades_nextstrain.tmp 2>&1 | tee {log}
        python scripts/overwrite_recombinant_clades.py \
            --clades clades_nextstrain.tmp \
            --internal-pango {input.internal_pango} \
            --alias {input.alias} \
            --clade-type clade_nextstrain \
            --output {output.node_data}
        rm clades_nextstrain.tmp
        cp {output.node_data} {output.node_data_nextstrain}
        sed -i'' 's/clade_membership/clade_nextstrain/gi' {output.node_data_nextstrain}
        """


rule clades_who:
    input:
        tree=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        nuc_muts=rules.ancestral.output.node_data,
        clades="builds/clades_who.tsv",
        internal_pango=rules.internal_pango.output.node_data,
        alias=rules.download_pango_alias.output,
    output:
        node_data=build_dir + "/{build_name}/clades_who.json",
    log:
        "logs/clades_{build_name}.txt",
    benchmark:
        "benchmarks/clades_{build_name}.txt"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data clades_who.tmp 2>&1 | tee {log}
        python scripts/overwrite_recombinant_clades.py \
            --clades clades_who.tmp \
            --internal-pango {input.internal_pango} \
            --alias {input.alias} \
            --clade-type clade_who \
            --output {output.node_data}
        rm clades_who.tmp
        sed -i'' 's/clade_membership/clade_who/gi' {output.node_data}
        """


rule download_designation_dates:
    output:
        designation_dates="builds/{build_name}/designation_dates.tsv",
    shell:
        """
        curl https://raw.githubusercontent.com/corneliusroemer/pango-designation-dates/main/data/lineage_designation_date.csv \
            | csv2tsv > {output}
        """


rule add_designation_date_to_meta:
    input:
        metadata="builds/{build_name}/metadata.tsv",
        designation_dates=rules.download_designation_dates.output.designation_dates,
    output:
        metadata="builds/{build_name}/metadata_with_designation_date.tsv",
    shell:
        """
        tsv-join -H --filter-file {input.designation_dates} \
            --key-fields 1 \
            --append-fields 2 \
            {input.metadata} \
        > {output}
        """


rule colors:
    input:
        ordering="defaults/color_ordering.tsv",
        color_schemes=config["files"]["color_schemes"],
        metadata="builds/{build_name}/metadata.tsv",
    output:
        colors=build_dir + "/{build_name}/colors.tsv",
    log:
        "logs/colors_{build_name}.txt",
    benchmark:
        "benchmarks/colors_{build_name}.txt"
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata} 2>&1 | tee {log}
        """


def _get_node_data_by_wildcards(wildcards):
    """Return a list of node data files to include for a given build's wildcards."""
    # Define inputs shared by all builds.
    wildcards_dict = dict(wildcards)
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.clades_legacy.output.node_data,
        rules.clades.output.node_data,
        rules.clades.output.node_data_nextstrain,
        rules.clades_who.output.node_data,
        rules.aa_muts_explicit.output.node_data,
        rules.internal_pango.output.node_data,
    ]

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]
    return inputs


rule export:
    input:
        tree=rules.refine.output.tree,
        metadata="builds/{build_name}/metadata_with_designation_date.tsv",
        node_data=_get_node_data_by_wildcards,
        auspice_config=lambda w: config["builds"][w.build_name]["auspice_config"]
        if "auspice_config" in config["builds"][w.build_name]
        else config["files"]["auspice_config"],
        description=lambda w: config["builds"][w.build_name]["description"]
        if "description" in config["builds"][w.build_name]
        else config["files"]["description"],
        colors=lambda w: rules.colors.output.colors.format(**w),
    output:
        auspice_json="auspice/{build_name}/auspice_raw.json",
        root_json="auspice/{build_name}/auspice_raw_root-sequence.json",
    log:
        "logs/export_{build_name}.txt",
    benchmark:
        "benchmarks/export_{build_name}.txt"
    params:
        title=lambda w: config["builds"][w.build_name].get(
            "title", "SARS-CoV-2 phylogeny"
        ),
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --title {params.title:q} \
            --description {input.description} \
            --include-root-sequence \
            --output {output.auspice_json} 2>&1 | tee {log};
        """


rule add_branch_labels:
    input:
        auspice_json=rules.export.output.auspice_json,
        mutations=rules.aa_muts_explicit.output.node_data,
    output:
        auspice_json="auspice/{build_name}/auspice.json",
    log:
        "logs/add_branch_labels_{build_name}.txt",
    shell:
        """
        python3 scripts/add_branch_labels.py \
            --input {input.auspice_json} \
            --mutations {input.mutations} \
            --output {output.auspice_json}
        """


rule remove_recombinants_from_auspice:
    input:
        auspice_json=rules.add_branch_labels.output.auspice_json,
    output:
        auspice_json="auspice/{build_name}/auspice_without_recombinants.json",
    log:
        "logs/remove_recombinants_from_auspice_{build_name}.txt",
    shell:
        """
        python3 scripts/remove_recombinants_from_auspice.py \
            --input {input.auspice_json} \
            --output {output.auspice_json}
        """


rule produce_trees:
    input:
        "auspice/nextclade/auspice.json",
        "auspice/nextclade/auspice_without_recombinants.json",
