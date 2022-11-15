import typer

# python scripts/overwrite_recombinant_clades.py \
# --clades {input.clades_json} \
# --output {output.clades_json}
def main(
    clades: str = "",
    internal_pango: str = "",
    alias: str = "",
    clade_type: str = "",
    output: str = "",
):
    import json
    from pango_aliasor.aliasor import Aliasor

    # Load clades.json
    with open(clades, "r") as f:
        clades = json.load(f)

    with open(internal_pango, "r") as f:
        internal_pango = json.load(f)
    
    aliasor = Aliasor(alias)

    # Overwrite values with `recombinant` where `key` starts with X
    # May want to set up explicit mapping for each Pango recombinant -> clade, WHO name etc
    for node, value in clades["nodes"].items():
        try:
            uncompressed = aliasor.uncompress(internal_pango["nodes"].get(node,{}).get("Nextclade_pango","")).split(".")[0]
        except:
            uncompressed = ""

        if uncompressed.startswith("X") or internal_pango["nodes"].get(node,{}).get("Nextclade_pango","").startswith(
            "X"
        ):
            if uncompressed not in ["XBB"]:
                value["clade_membership"] = "recombinant"
                value.pop("clade_annotation", None)
            elif clade_type == "clade_who":
                value["clade_membership"] = "Omicron"
        
        if clade_type != "clade_membership":
            value.pop("clade_annotation", None)

    # Write clades.json
    with open(output, "w") as f:
        json.dump(clades, f, indent=2)


if __name__ == "__main__":
    typer.run(main)
    # for debugging
    # main(clades="builds/nextclade/clades.json", output="builds/nextclade/clades_with_recombinants.json")
