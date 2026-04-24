#!/usr/bin/env python3
"""
Find ambiguous terminal zero-length branches in an Auspice v2 JSON tree.

This script treats branch length as the number of mutations recorded on the
child branch in `branch_attrs.mutations`. It reports internal nodes with two
or more terminal children whose branch length is zero.

Usage:
    python devel/check_zero_length_edges.py
    python devel/check_zero_length_edges.py auspice/BA.2.86/auspice.json
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any


DEFAULT_AUSPICE_JSON = Path("auspice/wuhan/auspice.json")


def attr_value(value: Any) -> Any:
    """Handle both raw Auspice attrs and {"value": ...} style attrs."""
    if isinstance(value, dict) and "value" in value:
        return value["value"]
    return value


def get_node_attr(node: dict[str, Any], key: str) -> Any:
    return attr_value(node.get("node_attrs", {}).get(key))


def get_label(node: dict[str, Any]) -> str:
    name = node.get("name", "?")
    pango = get_node_attr(node, "Nextclade_pango")
    clade = get_node_attr(node, "clade_membership")

    parts = [str(name)]
    if pango:
        parts.append(f"pango={pango}")
    if clade:
        parts.append(f"clade={clade}")
    return " ".join(parts)


def count_branch_mutations(node: dict[str, Any]) -> int:
    mutations = node.get("branch_attrs", {}).get("mutations", {})
    if not isinstance(mutations, dict):
        return 0
    return sum(len(values) for values in mutations.values() if isinstance(values, list))


def is_terminal(node: dict[str, Any]) -> bool:
    return not node.get("children")


def walk_nodes(root: dict[str, Any]):
    stack = [(root, 0)]
    while stack:
        node, depth = stack.pop()
        yield node, depth
        for child in reversed(node.get("children", [])):
            stack.append((child, depth + 1))


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Print internal nodes with multiple terminal children that have "
            "zero branch mutations in an Auspice JSON tree."
        )
    )
    parser.add_argument(
        "auspice_json",
        nargs="?",
        type=Path,
        default=DEFAULT_AUSPICE_JSON,
        help=f"Path to auspice.json (default: {DEFAULT_AUSPICE_JSON})",
    )
    parser.add_argument(
        "--fail-on-found",
        action="store_true",
        help="Exit with status 1 if any ambiguous zero-length terminal sibling groups are found.",
    )
    args = parser.parse_args()

    with args.auspice_json.open() as f:
        data = json.load(f)

    root = data.get("tree", data)
    internal_nodes = 0
    terminal_children = 0
    zero_terminal_children = 0
    violations: list[tuple[int, dict[str, Any], list[dict[str, Any]]]] = []

    for node, depth in walk_nodes(root):
        children = node.get("children", [])
        if not children:
            continue

        internal_nodes += 1
        zero_terminal_siblings = []
        for child in children:
            if not is_terminal(child):
                continue
            terminal_children += 1
            if count_branch_mutations(child) == 0:
                zero_terminal_children += 1
                zero_terminal_siblings.append(child)

        if len(zero_terminal_siblings) >= 2:
            violations.append((depth, node, zero_terminal_siblings))

    print(f"Input: {args.auspice_json}")
    print(f"Internal nodes checked: {internal_nodes}")
    print(f"Terminal children checked: {terminal_children}")
    print(f"Zero-length terminal children: {zero_terminal_children}")
    print(f"Internal nodes with >=2 zero-length terminal children: {len(violations)}")

    if not violations:
        return 0

    print()
    print("depth\tzero_terminal_children\tparent\tchildren")
    for depth, parent, children in violations:
        child_labels = " | ".join(get_label(child) for child in children)
        print(f"{depth}\t{len(children)}\t{get_label(parent)}\t{child_labels}")

    return 1 if args.fail_on_found else 0


if __name__ == "__main__":
    sys.exit(main())
