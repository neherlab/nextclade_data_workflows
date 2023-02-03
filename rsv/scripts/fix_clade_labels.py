import json, argparse

def replace_clade_recursive(node, label_key):
    if label_key in node["node_attrs"]:
        if "labels" not in node["branch_attrs"]:
            node["branch_attrs"]["labels"] = {}
        node["branch_attrs"]["labels"][label_key] = node["node_attrs"][label_key]["value"]
        node["node_attrs"].pop(label_key)
    if "children" in node:
        for child in node["children"]:
            replace_clade_recursive(child, label_key)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="fix genome clade info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-auspice-json', type=str, required=True, help="input auspice_json")
    parser.add_argument('--label-keys', type=str, nargs='+', required=True, help="input auspice_json")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()
    print(args.label_keys)
    with open(args.input_auspice_json, 'r') as fh:
        data = json.load(fh)

    for label_key in args.label_keys:
        data["meta"]["colorings"] = [x for x in data["meta"]["colorings"]
                                    if x["key"] != label_key]
        replace_clade_recursive(data['tree'], label_key)

    with open(args.output, 'w') as fh:
        json.dump(data, fh, indent=0)
