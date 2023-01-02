import click
import json

def swap(node):
    if "children" in node:
        for child in node["children"]:
            swap(child)
    elif "strainName" in node["node_attrs"]:
        accession = node["name"]
        node["name"] = node["node_attrs"]["strainName"]["value"]
        node["node_attrs"]["accessionID"] = {"value":  accession}
        del node["node_attrs"]["strainName"]
    else:
        print(f"Error in node:{node}")

def add_fake_clade_recursive(node, add_fake_clade):
    if "children" in node:
        for child in node["children"]:
            add_fake_clade_recursive(child, add_fake_clade)

    node["node_attrs"]["clade_membership"] = {"value": add_fake_clade}


@click.command()
@click.option("--tree-input-path", required=True)
@click.option("--tree-output-path", required=True)
@click.option("--add-fake-clade", required=False)
def main(tree_input_path, tree_output_path, add_fake_clade=None):
    with open(tree_input_path) as f:
        tree = json.load(f)

    swap(tree["tree"])
    if add_fake_clade:
        add_fake_clade_recursive(tree["tree"], add_fake_clade)

    for coloring in tree["meta"]["colorings"]:
        if coloring["key"] == "strainName":
            coloring["key"] = "accessionID"
            coloring["title"] = "Accession"

    with open(tree_output_path, "w") as f:
        json.dump(tree, f, indent=2)

if __name__ == "__main__":
    main()
