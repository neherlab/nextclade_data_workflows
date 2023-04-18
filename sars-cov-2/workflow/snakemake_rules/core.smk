localrules:
    add_branch_labels,
    colors,
    internal_pango,
    add_recombinants_to_tree,
    generate_nextclade_ba2_tsv,
    generate_nextclade_wuhan_tsv,
    download_nextclade_dataset,


genes = [
    "ORF1a",
    "ORF1b",
    "S",
    "ORF3a",
    "M",
    "N",
    "E",
    "ORF6",
    "ORF7a",
    "ORF7b",
    "ORF8",
    "ORF9b",
]


rule generate_nextclade_ba2_tsv:
    input:
        sequences="builds/{build_name}/sequences.fasta",
    output:
        tsv="builds/{build_name}/nextclade_ba2.tsv",
    shell:
        """
        # nextclade run -d sars-cov-2-21L -a auspice/21L/auspice.json {input.sequences} -t {output.tsv}
        nextclade run -d sars-cov-2-21L {input.sequences} -t {output.tsv}
        """


rule generate_nextclade_wuhan_tsv:
    input:
        sequences="builds/{build_name}/sequences.fasta",
    output:
        tsv="builds/{build_name}/nextclade_wuhan.tsv",
    shell:
        """
        nextclade run -d sars-cov-2 {input.sequences} -t {output.tsv}
        """


rule add_nextclade_columns_to_meta:
    input:
        metadata="builds/{build_name}/metadata_with_designation_date.tsv",
        ba2_tsv=rules.generate_nextclade_ba2_tsv.output.tsv,
        wuhan_tsv=rules.generate_nextclade_wuhan_tsv.output.tsv,
    output:
        metadata="builds/{build_name}/metadata_with_bloom_scores.tsv",
    shell:
        """
        tsv-join -H --filter-file {input.ba2_tsv} \
            --key-fields seqName \
            --data-fields strain \
            --append-fields immune_escape,ace2_binding \
            {input.metadata} \
        | tsv-join -H --filter-file {input.wuhan_tsv} \
            --key-fields seqName \
            --data-fields strain \
            --append-fields qc.overallScore,totalFrameShifts,qc.stopCodons.totalStopCodons,totalAminoacidSubstitutions,totalAminoacidDeletions \
        > {output}
        """


rule align:
    """
    Only used for translations
    As `sequences.fasta` is already aligned
    """
    input:
        sequences="builds/{build_name}/sequences.fasta",
        genemap="defaults/annotation.gff",
        reference="defaults/reference_seq.fasta",
    output:
        alignment="builds/{build_name}/aligned.fasta",
        translations=expand(
            "builds/{{build_name}}/translations/aligned.gene.{genes}.fasta",
            genes=genes,
        ),
    params:
        outdir=lambda w: f"builds/{w.build_name}/translations/aligned.gene.{{gene}}.fasta",
        genes=",".join(genes),
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
            --output-fasta {output.alignment}
        """


rule mask:
    input:
        alignment="builds/{build_name}/sequences.fasta",
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


rule prune_constraint_tree:
    """
    Prune constraint tree to only include sequences in the alignment
    """
    input:
        constraint_tree="defaults/constraint.nwk",
        strains=rules.get_strains.output.strains,
    output:
        constraint_tree="builds/{build_name}/constraint.nwk",
    shell:
        """
        python3 scripts/prune_constraint_tree.py \
            --constraint-tree {input.constraint_tree} \
            --strains {input.strains} \
            --output {output.constraint_tree}
        """


rule tree:
    input:
        alignment="builds/{build_name}/masked_without_recombinants.fasta",
        constraint_tree="builds/{build_name}/constraint.nwk",
        exclude_sites="defaults/exclude_sites.tsv",
    output:
        tree="builds/{build_name}/tree_raw.nwk",
    params:
        args=lambda w, input: f"-czb -g {input.constraint_tree} -ninit 1 -n 1",
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --exclude-sites {input.exclude_sites} \
            --tree-builder-args {params.args:q} \
            --output {output.tree} \
            --nthreads {threads} 2>&1
        """


rule recombinant_tree:
    input:
        alignment="builds/{build_name}/masked_recombinant_{recombinant}.fasta",
    output:
        tree="builds/{build_name}/tree_raw_recombinant_{recombinant}.nwk",
    params:
        constraint=lambda w: (
            f"-g profiles/clades/constraint_{w.recombinant}.nwk"
            if os.path.exists(f"profiles/clades/constraint_{w.recombinant}.nwk")
            else ""
        ),
    shell:
        """
        # Check if there are 3 or more sequences in the alignment
        if [ $(grep -c ">" {input.alignment}) -lt 3 ]; then
            python scripts/simple_tree.py \
                --alignment {input.alignment} \
                --tree {output.tree};
        else 
            augur tree \
                --alignment {input.alignment} \
                --tree-builder-args "-czb -ninit 1 -n 1 {params.constraint}" \
                --nthreads 1 \
                --output {output.tree};
        fi
        """


rule add_recombinants_to_tree:
    """
    Adding recombinants to root of raw tree
    """
    input:
        tree=rules.tree.output.tree,
        recombinants="builds/{build_name}/recombinants.txt",
        recombinant_trees=expand(
            "builds/{{build_name}}/tree_raw_recombinant_{recombinant}.nwk",
            recombinant=config["tree-recombinants"],
        ),
    output:
        tree="builds/{build_name}/tree_with_recombinants.nwk",
    params:
        root=lambda w: config["root"][w.build_name],
        joined_trees=lambda w, input: ",".join(input.recombinant_trees),
    shell:
        """
        python scripts/add_recombinants.py \
            --tree {input.tree} \
            --recombinants {input.recombinants} \
            --recombinant-trees {params.joined_trees} \
            --root {params.root} \
            --output {output.tree}
        """


rule refine:
    input:
        tree=rules.add_recombinants_to_tree.output.tree,
        alignment="builds/{build_name}/sequences.fasta",
        metadata="builds/{build_name}/metadata.tsv",
    output:
        tree="builds/{build_name}/tree.nwk",
        node_data="builds/{build_name}/branch_lengths.json",
    params:
        root=lambda w: config["root"][w.build_name],
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --keep-root \
            --divergence-unit mutations
        """


rule ancestral:
    input:
        tree=rules.refine.output.tree,
        alignment="builds/{build_name}/sequences.fasta",
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
        reference="defaults/reference_seq.gb",
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
        translations=rules.align.output.translations,
    output:
        node_data="builds/{build_name}/aa_muts_explicit.json",
        translations=expand(
            "builds/{{build_name}}/translations/aligned.gene.{gene}_withInternalNodes.fasta",
            gene=genes,
        ),
    params:
        genes=" ".join(genes),
    shell:
        """
        python3 scripts/explicit_translation.py \
            --tree {input.tree} \
            --translations {input.translations:q} \
            --genes {params.genes}\
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
            --build-name {wildcards.build_name} \
            --output {output.node_data} \
            --field-name Nextclade_pango
        """


rule preprocess_clades:
    input:
        clades="builds/clades{clade_type}.tsv",
        outgroup="profiles/clades/{build_name}/outgroup.tsv",
    output:
        clades="builds/{build_name}/clades{clade_type}.tsv",
    wildcard_constraints:
        clade_type=".*",  # Snakemake wildcard default is ".+" which doesn't match empty strings
    shell:
        """
        cp {input.clades} {output.clades};
        cat <(echo) {input.outgroup} >> {output.clades};
        if [ {wildcards.build_name} = 21L ]; then
            for clade in 19A 19B 20A 20B 20C 20D 20E 20F 20G 20H 20I \
                20J 21A 21B 21C 21D 21E 21F 21G 21H 21I 21J 21K 21M \
                Alpha Beta Gamma Delta Epsilon Eta Theta Iota Kappa Lambda Mu;
            do
                sed -i "/$clade/d" {output.clades};
            done
        fi
        """


rule clades_legacy:
    input:
        tree=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        nuc_muts=rules.ancestral.output.node_data,
        clades="builds/{build_name}/clades.tsv",
        internal_pango=rules.internal_pango.output.node_data,
        alias=rules.download_pango_alias.output,
    output:
        node_data="builds/{build_name}/clades_legacy.json",
    params:
        tmp="builds/{build_name}/clades_legacy.tmp",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {params.tmp}
        python scripts/overwrite_recombinant_clades.py \
            --clades {params.tmp} \
            --internal-pango {input.internal_pango} \
            --alias {input.alias} \
            --clade-type clade_legacy \
            --output {output.node_data}
        rm {params.tmp}
        sed -i'' 's/clade_membership/clade_legacy/gi' {output.node_data}
        """


rule clades:
    input:
        tree=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        nuc_muts=rules.ancestral.output.node_data,
        clades="builds/{build_name}/clades_nextstrain.tsv",
        internal_pango=rules.internal_pango.output.node_data,
        alias=rules.download_pango_alias.output,
    output:
        node_data="builds/{build_name}/clades.json",
        node_data_nextstrain="builds/{build_name}/clades_nextstrain.json",
    params:
        tmp="builds/{build_name}/clades_nextstrain.tmp",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {params.tmp}
        python scripts/overwrite_recombinant_clades.py \
            --clades {params.tmp} \
            --internal-pango {input.internal_pango} \
            --alias {input.alias} \
            --clade-type clade_nextstrain \
            --output {output.node_data}
        rm {params.tmp}
        cp {output.node_data} {output.node_data_nextstrain}
        sed -i'' 's/clade_membership/clade_nextstrain/gi' {output.node_data_nextstrain}
        """


rule clades_who:
    input:
        tree=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        nuc_muts=rules.ancestral.output.node_data,
        clades="builds/{build_name}/clades_who.tsv",
        internal_pango=rules.internal_pango.output.node_data,
        alias=rules.download_pango_alias.output,
    output:
        node_data="builds/{build_name}/clades_who.json",
    params:
        tmp="builds/{build_name}/clades_who.tmp",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {params.tmp}
        python scripts/overwrite_recombinant_clades.py \
            --clades {params.tmp} \
            --internal-pango {input.internal_pango} \
            --alias {input.alias} \
            --clade-type clade_who \
            --output {output.node_data}
        rm {params.tmp}
        sed -i'' 's/clade_membership/clade_who/gi' {output.node_data}
        sed -i'' 's/unassigned//gi' {output.node_data} # replace unassigned with empty string
        """


rule colors:
    input:
        ordering="defaults/color_ordering.tsv",
        color_schemes="defaults/color_schemes.tsv",
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
        metadata=rules.add_nextclade_columns_to_meta.output.metadata,
        node_data=_get_node_data_by_wildcards,
        colors=rules.colors.output.colors,
        auspice_config="profiles/clades/{build_name}/auspice_config.json",
        description="profiles/clades/description.md",
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
            --minify-json \
            --output {output.auspice_json}
        """


rule download_nextclade_dataset:
    """
    Download Nextclade dataset so Nextclade can run without internet connection
    """
    output:
        dataset="builds/{build_name}/nextclade_dataset.zip",
    params:
        dataset=lambda w: "sars-cov-2" if w.build_name == "wuhan" else "sars-cov-2-21L",
    shell:
        """
        nextclade dataset get \
            --name {params.dataset} \
            --output-zip {output.dataset}
        """


rule generate_priors:
    """
    Run nextclade on generated tree to get placement priors
    """
    input:
        fasta=rules.download_sequences.output,
        tree=rules.export.output.auspice_json,
        dataset=rules.download_nextclade_dataset.output.dataset,
    output:
        ndjson="builds/{build_name}/nextclade.ndjson.zst",
    shell:
        """
        zstdcat -T2 {input.fasta} | \
        seqkit sample -p 0.1 -w0 | \
        nextclade run \
            -D {input.dataset} \
            -a {input.tree} \
            --include-nearest-node-info \
            --output-ndjson {output.ndjson}
        """


rule get_nearest_nodes_from_ndjson:
    """
    Extract nearestNodes from Nextclade output
    """
    input:
        ndjson=rules.generate_priors.output.ndjson,
    output:
        ndjson="builds/{build_name}/nearest_nodes.ndjson",
    shell:
        """
        zstdcat {input.ndjson} | \
        jq -c '{{seqName,nearestNodes}}' > {output.ndjson}
        """


rule add_priors:
    """
    Update placement priors
    """
    input:
        auspice_json=rules.export.output.auspice_json,
        ndjson=rules.get_nearest_nodes_from_ndjson.output.ndjson,
    output:
        auspice_json="builds/{build_name}/auspice_priors.json",
    shell:
        """
        python3 scripts/add_priors.py \
            --tree {input.auspice_json} \
            --ndjson {input.ndjson} \
            --output {output.auspice_json}
        """


rule add_branch_labels:
    input:
        auspice_json=rules.add_priors.output.auspice_json,
        mutations=rules.aa_muts_explicit.output.node_data,
    output:
        auspice_json="auspice/{build_name}/auspice_max.json",
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
        auspice_json="auspice/{build_name}/auspice_without_recombinants_max.json",
    shell:
        """
        python3 scripts/remove_recombinants_from_auspice.py \
            --input {input.auspice_json} \
            --output {output.auspice_json}
        """


rule minify_json:
    input:
        "auspice/{build_name}/{build_type}_max.json",
    output:
        "auspice/{build_name}/{build_type}.json",
    shell:
        """
        jq -c . {input} > {output}
        """


rule produce_trees:
    input:
        "auspice/wuhan/auspice.json",
        "auspice/wuhan/auspice_without_recombinants.json",
        "auspice/21L/auspice.json",
        "auspice/21L/auspice_without_recombinants.json",
