# Nextclade reference tree workflow for SARS-CoV-2

Builds the Nextclade SARS-CoV-2 reference trees and dataset files for the three
datasets published in `nextclade_data`:

| build      | dataset path in `nextclade_data`                 | staging tree on S3                   |
| ---------- | ------------------------------------------------ | ------------------------------------ |
| `wuhan`    | `data/nextstrain/sars-cov-2/wuhan-hu-1/orfs`     | `nextclade_sars-cov-2.json`          |
| `proteins` | `data/nextstrain/sars-cov-2/wuhan-hu-1/proteins` | `nextclade_sars-cov-2_proteins.json` |
| `BA.2.86`  | `data/nextstrain/sars-cov-2/BA.2.86`             | `nextclade_sars-cov-2_BA.2.86.json`  |

The `XBB` and `BA.2` datasets in `nextclade_data` are frozen — they are not built
here and not updated on release.

## Running the workflow

```sh
just build   # snakemake -c10 --profile profiles/clades -F
just retry   # snakemake -c10 --profile profiles/clades --ri
```

Use the `py13` environment (`micromamba run -n py13 snakemake ...`) when invoking
snakemake directly, and `-R <rule>` to force a rerun of a single rule.

The `all` rule deploys each build's tree to `s3://nextstrain-staging` and writes
the assembled dataset files to `datasets/{path}/`.

## Sanity checks after a build

```sh
uv run devel/check_tree_lineages.py datasets/wuhan-hu-1/orfs/tree.json   # parent/child lineage consistency
python devel/check_zero_length_edges.py auspice/wuhan/auspice.json       # ambiguous zero-length terminals
```

## Maintainer checklist: adding new clades

When `nextstrain/ncov` designates a new Nextstrain clade, every one of these has
to be touched. Grepping for the previous JN.1-descendant clade name (e.g. `26F`)
finds them all: outside generated files and `defaults/labeled_muts.json`, only
`color_ordering.tsv`, `generate_virus_properties.py` and the three
`auspice_config.json` files name clades.

1. **Clade definitions** come from `nextstrain/ncov` (`defaults/clades.tsv` plus
   clade display names) via the preprocess rules that generate
   `builds/clades_nextstrain.tsv` and `builds/clades_display.tsv`. They are
   fetched from the ncov branch named by `ncovbranch` in
   `profiles/clades/builds.yaml`; point it at the designation branch while the
   ncov PR is unmerged and set it back to `master` once merged.
2. **Recombinant allow-list** — if the clade is recombinant, add its unaliased
   prefix to the list in `scripts/overwrite_recombinant_clades.py`
   (currently `XBB`, `XDV.1`, `XEC`, `XFG`, `XFC`, `XFJ`), otherwise the clade is
   collapsed to `recombinant`. A sublineage clade is already covered by its
   parent's prefix — 26A (XFG.1.1) needs nothing beyond `XFG`.
3. **`defaults/color_ordering.tsv`** — add the clade to *all four* clade
   sections: `clade_membership`, `clade`, `clade_legacy` (as `26A (Omicron)`) and
   `clade_nextstrain`. The colors themselves are assigned positionally by
   `scripts/assign-colors.py` from `defaults/color_schemes.tsv`, so there are no
   hex codes to pick; adding entries shifts the palette for all clades.
4. **`profiles/clades/*/auspice_config.json`** (all three) — the `ref_nodes`
   search criteria (JN.1, XBB.1.5) restrict each reference node to the query
   clades descending from it. Add a new clade to the `qry` clade list of the
   reference node it descends from (in practice JN.1).
5. **Labeled mutations** — see below. A new clade has no characteristic
   mutations until it appears in the ncov metadata, which only happens *after*
   this dataset is released, so the first pass needs a manual override.
6. **Verify after the build**: the clade appears in `builds/*/clades_nextstrain.json`,
   `builds/*/colors.tsv` and `auspice/*/auspice.json`, and is not collapsed to
   `recombinant`.

### Known gap

`defaults/sequences.fasta` — the example sequences copied verbatim into every
published dataset — has not been refreshed since 2023-12-05. All 165 sequences
predate every clade designated since; they should be resampled to include
currently circulating clades.

## Labeled mutations / virus properties

`scripts/generate_virus_properties.py` computes characteristic mutations per
clade from the full ncov metadata. It needs the private metadata dump in the
working directory first (~1.6 GB; the run itself takes several minutes):

```sh
aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.zst .
micromamba run -n py12 python scripts/generate_virus_properties.py
```

It writes two files with the same mutation map:

- `defaults/labeled_muts.json` — consumed by the build via
  `scripts/assemble_pathogen_json.py` and baked into `builds/*/pathogen.json`
- `virus_properties.json` — the same map plus a reverse index, used to patch
  already-released `pathogen.json` files

The metadata's `clade_nextstrain` column is produced by a Nextclade run against
the *released* dataset, so a clade designated since the last release is absent
from it. Until then, override the assignment near the top of the script by
unaliased lineage, e.g.

```python
df.loc[is_lineage_or_sublineage(df["unaliased"], "XFG.1.1"), "Nextstrain_clade"] = "26A"
```

Match exact-or-sublineage, not a bare `str.startswith`: `startswith("XFG.1.1")`
would also match `XFG.1.10`.

To patch released datasets without a rebuild, use
`scripts/update_virus_properties.py` (merges `virus_properties.json` into the
`mutLabels` of each `nextclade_data` `pathogen.json`), or equivalently:

```sh
cd ~/code
for d in BA.2.86 wuhan-hu-1/{orfs,proteins}; do
  F=./nextclade_data/data/nextstrain/sars-cov-2/$d/pathogen.json
  jq --slurpfile v ./nextclade_data_workflows/sars-cov-2/virus_properties.json '.mutLabels.nucMutLabelMap=$v[0].nucMutLabelMap' $F | sponge $F
done
```

## Cutting a new dataset release

1. **Update the changelog.** Pull `~/code/pango-designation-dates` first, then
   regenerate the `## Unreleased` section of `defaults/CHANGELOG.md`, starting
   the day after the previous release's window ended and ending the day before
   today (designations made today are not in the tree yet):

   ```sh
   python devel/update_changelog_lineages.py 2026-06-09 2026-09-06
   ```

   Add a line by hand for any newly designated Nextstrain clades, linking the
   ncov PRs, e.g.

   > - Include the 2 newly designated Nextstrain clades 26A (XFG.1.1) and 26B
   >   (XFG.3.4.1). See <https://github.com/nextstrain/ncov/pull/1204> and
   >   <https://github.com/nextstrain/ncov/pull/1206> for details.

2. **Check the staging trees** are the build you want
   (`aws s3 ls s3://nextstrain-staging/`; the tree's `meta.updated` carries the
   build date).

3. **Copy trees + changelog into `nextclade_data`** and publish to prod:

   ```sh
   python devel/release.py
   ```

   This downloads each staging tree to
   `~/code/nextclade_data/data/nextstrain/sars-cov-2/<path>/tree.json`, prepends
   `defaults/CHANGELOG.md` to each dataset's `CHANGELOG.md`, copies
   `builds/<build>/pathogen.json`, and copies the staging tree to
   `s3://nextstrain-data`. Skip the `pathogen.json` copy when nothing that feeds
   it (`defaults/labeled_muts.json`, `defaults/pathogen_json_base.json`,
   `profiles/clades/*/pathogen_json/extras.json`) has changed since the last
   release.

   Leave the `## Unreleased` heading in place — `nextclade_data`'s release CI
   replaces it with the release timestamp.

4. **Open a PR** in `nextclade_data` against `master` with the three `tree.json`
   and three `CHANGELOG.md` files.

## Immune escape and ACE2 phenotype data

The datasets no longer include RBD immune escape or ACE2 affinity scores:
`profiles/clades/*/pathogen_json/{escape,ace2}.json` are empty (`{}`), so
`scripts/assemble_pathogen_json.py` adds no `phenotypeData`, and no
`auspice_config.json` defines the `immune_escape`/`ace2_binding` colorings.

## TODO

- Refresh `defaults/sequences.fasta` (see "Known gap" above) (create script to add for clades that don't yet have seqs)
