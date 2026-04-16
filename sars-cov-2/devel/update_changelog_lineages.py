#!/usr/bin/env python3
"""
Update the Unreleased section of defaults/CHANGELOG.md with newly designated
Pango lineages.

Usage:
    python update_changelog_lineages.py <start_date> [<end_date>]

    start_date: ISO date (YYYY-MM-DD), include lineages designated on or after this date
    end_date:   ISO date (YYYY-MM-DD), include lineages designated on or before this date
                (default: today)
"""

import csv
import re
import sys
from datetime import date, datetime
from pathlib import Path

CSV_PATH = Path("/Users/cr/code/pango-designation-dates/data/lineage_designation_date.csv")
CHANGELOG_PATH = Path(__file__).parent / "defaults" / "CHANGELOG.md"


def parse_date(s: str) -> date:
    return datetime.strptime(s, "%Y-%m-%d").date()


def load_lineages(start: date, end: date) -> list[tuple[str, date]]:
    results = []
    with CSV_PATH.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw = row["designation_date"].strip()
            if not raw:
                continue
            d = parse_date(raw)
            if start <= d <= end:
                results.append((row["lineage"].strip(), d))
    results.sort(key=lambda x: x[1])  # stable sort: preserve CSV order within each date
    return results


def build_list(lineages: list[tuple[str, date]]) -> str:
    return "\n".join(f"- {name} ({d})" for name, d in lineages)


def update_changelog(start: date, end: date, lineages: list[tuple[str, date]]) -> None:
    text = CHANGELOG_PATH.read_text()
    count = len(lineages)
    lineage_list = build_list(lineages)

    # Replace the summary line (count + date range); count may be a number or placeholder like XX
    text = re.sub(
        r"- Add all [\dXx]+ Pango lineages newly designated between .+? and .+?\.",
        f"- Add all {count} Pango lineages newly designated between {start} and {end}.",
        text,
    )

    # Replace the lineage list inside the <details> block
    text = re.sub(
        r"(<details>.*?<summary>[^<]*</summary>\s*\n)(.*?)(\n</details>)",
        lambda m: m.group(1) + lineage_list + "\n" + m.group(3),
        text,
        flags=re.DOTALL,
    )

    CHANGELOG_PATH.write_text(text)
    print(f"Updated {CHANGELOG_PATH} with {count} lineages ({start} to {end})")


def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    start = parse_date(sys.argv[1])
    end = parse_date(sys.argv[2]) if len(sys.argv) > 2 else date.today()

    lineages = load_lineages(start, end)
    if not lineages:
        print(f"No lineages found between {start} and {end}", file=sys.stderr)
        sys.exit(1)

    update_changelog(start, end, lineages)


if __name__ == "__main__":
    main()
