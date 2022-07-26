import typer

def main(input: str = "", output: str = ""):
    import json
    tree = json.load(open(input))
    tree["tree"]["children"][:] = [
        child for child in tree["tree"]["children"] if child["name"] != "rec_parent"
    ]
    with open(output, 'w') as f:
        json.dump(tree, f, indent=2)

if __name__ == '__main__':
    typer.run(main)
