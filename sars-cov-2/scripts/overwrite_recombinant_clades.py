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

    def startswith_in_list(string, startswith_list):
        for startswith in startswith_list:
            if string.startswith(startswith):
                return True
    for node, value in clades["nodes"].items():
        pango_lineage = internal_pango["nodes"].get(node,{}).get("Nextclade_pango","")
        try:
            unaliased = aliasor.uncompress(pango_lineage)
        except:
            unaliased = ""

        if unaliased.startswith("X"):
            # List all recombinants that define clade here
            # Needs to be updated if new recombinants become clades
            # XBB, XDV.1 are the only ones so far
            if not startswith_in_list(unaliased, ["XBB", "XDV.1", "XEC"]):
                value["clade_membership"] = "recombinant"
                value.pop("clade_annotation", None)
        
        if clade_type != "clade_nextstrain": # use clade_nextstrain for branch labels
            value.pop("clade_annotation", None)

    # Write clades.json
    with open(output, "w") as f:
        json.dump(clades, f, indent=2)


if __name__ == "__main__":
    typer.run(main)
    # for debugging
    # main(clades="builds/nextclade/clades.json", output="builds/nextclade/clades_with_recombinants.json")
