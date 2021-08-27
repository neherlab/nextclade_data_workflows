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

@click.command()
@click.option("--tree-input-path", required=True)
@click.option("--tree-output-path", required=True)
def main(tree_input_path, tree_output_path):
    with open(tree_input_path) as f:
        tree = json.load(f)

    swap(tree["tree"])

    for coloring in tree["meta"]["colorings"]:
        if coloring["key"] == "strainName":
            coloring["key"] = "accessionID"
            coloring["title"] = "Accession"
    
    with open(tree_output_path, "w") as f:
        json.dump(tree, f, indent=2)

if __name__ == "__main__":
    main()
