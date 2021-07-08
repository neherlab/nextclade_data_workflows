How to assign to strain: Yam/Vic
Mutation summary: look for nearest neighbour

Todo: extend to other segments

Refine: What should be run

QC: How to get rid of far away 


Sanitize script: tree and cut off long terminals
Use tree time: give it tree.

How to filter long branches?
Internal branch: little more conservative
Phylotree objects
Parameter: clock filter

Sample: monat region

Are insertions basically thrown away?

Nextalign, explicit translation

Why is date not shown on auspice?

Strain name vs accession number

What to subsample on? Year? Continent?

Are there clades I should use?

Download all data by year and join together

Color by metadata: Includes metadata in auspice.json



- Download right sequences
- Parse metadata
- Subsample
- Build tree

Make a tree for all virus types A, B, C and D


Allow specification for various settings of recency

Flu base from nextstrain: https://github.com/nextstrain/seasonal-flu/blob/master/Snakefile_base

Group by doesn't really equalise differences
Confused about sampling
Ambiguous dates excluded by default with min date?
Cannot specify which column to use?
Pandas filter does not work