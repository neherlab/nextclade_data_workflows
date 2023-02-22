localrules:
    add_branch_labels,
    colors,
    internal_pango,
    overwrite_recombinant_clades,
    add_recombinants_to_tree,
    remove_recombinants_from_alignment,
    identify_recombinants,


rule align:
    """
    TODO: Unnecessary, because synthetic sequences are already aligned
    """
    input:
        sequences="builds/{build_name}/sequences.fasta",
        genemap="defaults/annotation.gff",
        reference="defaults/reference_seq.fasta",
    output:
        alignment="builds/{build_name}/aligned.fasta",
        translations="builds/{build_name}/translations/aligned.gene.S.fasta",
    params:
        outdir=lambda w: f"builds/{w.build_name}/translations/aligned.gene.{{gene}}.fasta",
    threads: 4
    shell:
        """
        nextalign run \
            --jobs={threads} \
            --input-ref {input.reference} \
            --input-gene-map {input.genemap} \
            {input.sequences} \
            --output-translations {params.outdir} \
            --output-fasta {output.alignment} \
            2>&1
        """


rule mask:
    input:
        alignment=rules.align.output.alignment,
    output:
        alignment="builds/{build_name}/masked.fasta",
    shell:
        """
        augur mask \
            --sequences {input.alignment} \
            --mask-from-beginning 100 \
            --mask-from-end 100 \
            --mask-invalid \
            --output {output.alignment}
        """


rule separate_recombinants:
    input:
        alignment=rules.mask.output.alignment,
        alias_json=rules.download_pango_alias.output,
    output:
        without_recombinants="builds/{build_name}/masked_without_recombinants.fasta",
        recombinant_alignments=expand(
            "builds/{{build_name}}/masked_recombinant_{tree_recombinants}.fasta",
            tree_recombinants=config["tree-recombinants"],
        ),
        recombinants="builds/{build_name}/recombinants.txt",
    params:
        tree_recombinants=",".join(config["tree-recombinants"]),
    shell:
        """
        python3 scripts/separate_recombinants.py \
            --alignment {input.alignment} \
            --output-without-recombinants {output.without_recombinants} \
            --output-recombinant "builds/{wildcards.build_name}" \
            --alias-json {input.alias_json} \
            --tree-recombinants {params.tree_recombinants} \
            --recombinants {output.recombinants}
        """


rule tree:
    input:
        alignment="builds/nextclade/masked_without_recombinants.fasta",
        constraint_tree="defaults/constraint.nwk",
        exclude_sites="defaults/exclude_sites.tsv",
    output:
        tree="builds/nextclade/tree_raw.nwk",
    params:
        args=lambda w, input: f"'-czb -g {input.constraint_tree}'",
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
        tree="builds/{build_name}/tree.nwk",
        node_data="builds/{build_name}/branch_lengths.json",
    params:
        root=config["refine"]["root"],
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --keep-root \
            --divergence-unit mutations | tee {log}
        """


rule ancestral:
    input:
        tree=rules.refine.output.tree,
        alignment=rules.align.output.alignment,
    output:
        node_data="builds/{build_name}/nt_muts.json",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference joint \
            --infer-ambiguous
        """


rule translate:
    input:
        tree=rules.refine.output.tree,
        node_data=rules.ancestral.output.node_data,
        reference=config["files"]["reference"],
    output:
        node_data="builds/{build_name}/aa_muts.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data}
        """


rule aa_muts_explicit:
    input:
        tree=rules.refine.output.tree,
        translations=lambda w: rules.align.output.translations,
    output:
        node_data="builds/{build_name}/aa_muts_explicit.json",
        translations="builds/{build_name}/translations/aligned.gene.S_withInternalNodes.fasta",
    shell:
        """
        python3 scripts/explicit_translation.py \
            --tree {input.tree} \
            --translations {input.translations:q} \
            --genes S \
            --output {output.node_data}
        """


rule internal_pango:
    input:
        tree=rules.refine.output.tree,
        alias=rules.download_pango_alias.output,
        synthetic=rules.synthetic_pick.output,
        designations=rules.pango_strain_rename.output.pango_designations,
    output:
        node_data="builds/{build_name}/internal_pango.json",
    shell:
        """
        python scripts/internal_pango.py \
            --tree {input.tree} \
            --synthetic {input.synthetic} \
            --alias {input.alias} \
            --designations {input.designations} \
            --output {output.node_data} \
            --field-name Nextclade_pango
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
        node_data="builds/{build_name}/clades_legacy.json",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data clades_raw.tmp
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
        node_data="builds/{build_name}/clades.json",
        node_data_nextstrain="builds/{build_name}/clades_nextstrain.json",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data clades_nextstrain.tmp
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
        node_data="builds/{build_name}/clades_who.json",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data clades_who.tmp
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
        colors="builds/{build_name}/colors.tsv",
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata}
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
        auspice_config="profiles/clades/auspice_config.json",
        description="profiles/clades/description.md",
        colors=rules.colors.output.colors,
    output:
        auspice_json="auspice/{build_name}/auspice_raw.json",
        root_json="auspice/{build_name}/auspice_raw_root-sequence.json",
    params:
        title="SARS-CoV-2 phylogeny",
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
