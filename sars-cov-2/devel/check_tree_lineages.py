#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "pango-aliasor",
# ]
# ///
"""
Check tree.json for parent/child lineage consistency.
For each parent/child pair, fully unaliases lineages using pango_aliasor
and verifies that the parent's unaliased lineage is a logical prefix of
the child's unaliased lineage.

Usage:
    uv run check_tree_lineages.py datasets/wuhan-hu-1/orfs/tree.json
"""

import json
import sys
import argparse
from pango_aliasor.aliasor import Aliasor


def get_pango(node):
    attrs = node.get("node_attrs", {})
    pa = attrs.get("Nextclade_pango", {})
    if isinstance(pa, dict):
        return pa.get("value", "") or ""
    return str(pa) if pa else ""


def is_logical_child(parent_lineage, child_lineage):
    """
    Check if child_lineage is a logical descendant of parent_lineage.
    Both should be fully unaliased (e.g. B.1.1.529.2.86.1.1).
    """
    if not parent_lineage or not child_lineage:
        return None  # can't determine
    if parent_lineage == child_lineage:
        return True  # same lineage label (internal node)
    if child_lineage.startswith(parent_lineage + "."):
        return True
    return False


def levels_skipped(parent_lineage, child_lineage):
    """
    If child is a descendant of parent, return how many levels are skipped
    (0 = direct child, 1 = grandchild skipping one level, etc.).
    Returns None if not a descendant or same lineage.
    """
    if not parent_lineage or not child_lineage:
        return None
    if parent_lineage == child_lineage:
        return None
    if not child_lineage.startswith(parent_lineage + "."):
        return None
    parent_parts = parent_lineage.split(".")
    child_parts = child_lineage.split(".")
    return len(child_parts) - len(parent_parts) - 1


def lineage_depth(lineage):
    """Number of dot-separated components (e.g. B=1, B.1=2, B.1.1=3)."""
    if not lineage:
        return 0
    return len(lineage.split("."))


def walk_tree(node, aliasor, parent=None, parent_unaliased=None, violations=None, skipped_level=None, stats=None):
    if violations is None:
        violations = []
    if skipped_level is None:
        skipped_level = []
    if stats is None:
        stats = {"checked": 0, "skipped": 0, "violations": 0, "skipped_level": 0}

    pango = get_pango(node)
    node_unaliased = aliasor.uncompress(pango) if pango else ""
    node_name = node.get("name", "?")

    if parent is not None:
        parent_name = parent.get("name", "?")
        result = is_logical_child(parent_unaliased, node_unaliased)
        if result is None:
            stats["skipped"] += 1
        elif result is False:
            stats["checked"] += 1
            stats["violations"] += 1
            parent_pango = get_pango(parent)
            violations.append({
                "parent_name": parent_name,
                "parent_pango": parent_pango,
                "parent_unaliased": parent_unaliased,
                "parent_depth": lineage_depth(parent_unaliased),
                "child_name": node_name,
                "child_pango": pango,
                "child_unaliased": node_unaliased,
            })
        else:
            stats["checked"] += 1
            skipped = levels_skipped(parent_unaliased, node_unaliased)
            if skipped is not None and skipped > 0:
                stats["skipped_level"] += 1
                parent_pango = get_pango(parent)
                skipped_level.append({
                    "parent_name": parent_name,
                    "parent_pango": parent_pango,
                    "parent_unaliased": parent_unaliased,
                    "parent_depth": lineage_depth(parent_unaliased),
                    "child_name": node_name,
                    "child_pango": pango,
                    "child_unaliased": node_unaliased,
                    "levels_skipped": skipped,
                })

    for child in node.get("children", []):
        walk_tree(child, aliasor, parent=node, parent_unaliased=node_unaliased,
                  violations=violations, skipped_level=skipped_level, stats=stats)

    return violations, skipped_level, stats


def main():
    parser = argparse.ArgumentParser(description="Check tree.json parent/child lineage consistency")
    parser.add_argument("tree_json", help="Path to tree.json file")
    parser.add_argument("--max", type=int, default=100, help="Max violations to show (default: 100)")
    parser.add_argument("--all", action="store_true", help="Show all violations")
    parser.add_argument("--min-parent-depth", type=int, default=2, metavar="N",
                        help="Skip violations where parent lineage has fewer than N components "
                             "(e.g. 2 skips single-letter parents like B, A; default: 2)")
    args = parser.parse_args()

    aliasor = Aliasor()

    print(f"Loading {args.tree_json}...", file=sys.stderr)
    with open(args.tree_json) as f:
        data = json.load(f)

    tree = data.get("tree", data)
    print("Walking tree...", file=sys.stderr)
    violations, skipped_level, stats = walk_tree(tree, aliasor)

    filtered = [v for v in violations if v["parent_depth"] >= args.min_parent_depth]
    filtered_sl = [v for v in skipped_level if v["parent_depth"] >= args.min_parent_depth]
    skipped_shallow = len(violations) - len(filtered)

    print(f"Pairs checked: {stats['checked']}, skipped (no lineage): {stats['skipped']}, "
          f"wrong-parent violations: {stats['violations']}, "
          f"skipped-level violations: {stats['skipped_level']} "
          f"({skipped_shallow} filtered out with parent depth < {args.min_parent_depth})\n")

    DIM = "\033[2m"
    BOLD = "\033[1m"
    RESET = "\033[0m"

    def highlight_divergence(a, b):
        """Dim the common leading components, bold the diverging tail."""
        a_parts = a.split(".")
        b_parts = b.split(".")
        common = 0
        for x, y in zip(a_parts, b_parts):
            if x == y:
                common += 1
            else:
                break
        def fmt(parts):
            prefix = ".".join(parts[:common])
            suffix = ".".join(parts[common:])
            if prefix and suffix:
                return f"{DIM}{prefix}{RESET}.{BOLD}{suffix}{RESET}"
            elif prefix:
                return f"{DIM}{prefix}{RESET}"
            else:
                return f"{BOLD}{suffix}{RESET}"
        return fmt(a_parts), fmt(b_parts)

    c = [15, 15]
    header = f"{'Parent':<{c[0]}} {'Child':<{c[1]}}  Parent unaliased / Child unaliased"

    # --- Category 1: wrong parent ---
    if not filtered:
        print("No wrong-parent violations found after filtering.")
    else:
        limit = len(filtered) if args.all else min(args.max, len(filtered))
        print(f"=== Wrong-parent violations: {limit} of {len(filtered)} shown ===\n")
        print(header)
        print("-" * 80)
        for v in filtered[:limit]:
            pu, cu = highlight_divergence(v["parent_unaliased"], v["child_unaliased"])
            print(f"{v['parent_pango']:<{c[0]}} {v['child_pango']:<{c[1]}}  {pu}  /  {cu}")
        if limit < len(filtered):
            print(f"\n... and {len(filtered) - limit} more. Use --all or --max N.")

    print()

    # --- Category 2: skipped levels ---
    if not filtered_sl:
        print("No skipped-level violations found after filtering.")
    else:
        limit_sl = len(filtered_sl) if args.all else min(args.max, len(filtered_sl))
        print(f"=== Skipped-level violations (child not a direct child): {limit_sl} of {len(filtered_sl)} shown ===\n")
        print(header)
        print("-" * 80)
        for v in sorted(filtered_sl, key=lambda x: -x["levels_skipped"])[:limit_sl]:
            pu, cu = highlight_divergence(v["parent_unaliased"], v["child_unaliased"])
            skip_label = f"(+{v['levels_skipped']})"
            print(f"{v['parent_pango']:<{c[0]}} {v['child_pango']:<{c[1]}}  {pu}  /  {cu}  {skip_label}")
        if limit_sl < len(filtered_sl):
            print(f"\n... and {len(filtered_sl) - limit_sl} more. Use --all or --max N.")


if __name__ == "__main__":
    main()
