# Usage: python scripts/fix_rec_parent.py --tree tree.nwk --node-data muts.json --output muts_fixed.json [--node rec_parent]
# Pushes the mutations that augur ancestral put on the dummy rec_parent branch down onto its children
# (composing per site, nuc and every CDS), so rec_parent carries the root genotype and no node's genotype changes.
import argparse, json, re
from Bio import Phylo

MUT = re.compile(r"^(.)(\d+)(.)$")

def compose(upper, lower):
    st = {}
    for m in upper:
        a, p, b = MUT.match(m).groups(); st[int(p)] = [a, b]
    for m in lower:
        a, p, b = MUT.match(m).groups(); p = int(p)
        if p in st:
            if st[p][1] != a:
                raise ValueError(f"inconsistent chain at {p}: {st[p]} then {m}")
            st[p][1] = b
        else:
            st[p] = [a, b]
    return [f"{a}{p}{b}" for p, (a, b) in sorted(st.items()) if a != b]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True); ap.add_argument("--node-data", required=True)
    ap.add_argument("--output", required=True); ap.add_argument("--node", default="rec_parent")
    args = ap.parse_args()
    t = Phylo.read(args.tree, "newick")
    rp = next(c for c in t.find_clades() if c.name == args.node)
    parent = next(c for c in t.find_clades() if rp in c.clades)
    d = json.load(open(args.node_data)); nodes = d["nodes"]; rpd = nodes[args.node]
    n_moved = len(rpd.get("muts", [])) + sum(len(v) for v in rpd.get("aa_muts", {}).values())
    for c in rp.clades:
        cd = nodes[c.name]
        cd["muts"] = compose(rpd.get("muts", []), cd.get("muts", []))
        for g, ms in rpd.get("aa_muts", {}).items():
            cd.setdefault("aa_muts", {})[g] = compose(ms, cd.get("aa_muts", {}).get(g, []))
    rpd["muts"] = []
    rpd["aa_muts"] = {g: [] for g in rpd.get("aa_muts", {})}
    if "sequence" in rpd and "sequence" in nodes.get(parent.name, {}):
        rpd["sequence"] = nodes[parent.name]["sequence"]
    json.dump(d, open(args.output, "w"), indent=1)
    print(f"moved {n_moved} {args.node} mutations onto {len(rp.clades)} children")

if __name__ == "__main__":
    main()
