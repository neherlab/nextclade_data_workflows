# Nextclade reference tree workflow for SARS-CoV-2

## Running the workflow

```sh
snakemake -c10 --profile profiles/clades -F
```

## TODO

- Produce entire datasets in workflow, not just tree
  - pathogen.json needs labeled mutations, escape config
  - rest should be easy

## Creating a new dataset version

```sh
#!/bin/bash
OLD="2023-12-03T12:00:00Z"
NEW="2024-01-15T12:00:00Z"

SC2="sars-cov-2/references/MN908947"
SC221L="sars-cov-2-21L/references/BA.2"

declare -A EXTENSION=( ["$SC2"]="" ["$SC221L"]="_21L" )

for d in $SC2 $SC221L; do
  rm -rf ~/code/nextclade_data/data/datasets/$d/versions/$NEW
  cp -pr ~/code/nextclade_data/data/datasets/$d/versions/$OLD ~/code/nextclade_data/data/datasets/$d/versions/$NEW
  sed -i'.bak' "s/$OLD/$NEW/g" ~/code/nextclade_data/data/datasets/$d/versions/$NEW/files/tag.json;
  rm ~/code/nextclade_data/data/datasets/$d/versions/$NEW/files/tag.json.bak;
  aws s3 cp s3://nextstrain-staging/nextclade_sars-cov-2${EXTENSION[$d]}.json - | gzcat >~/code/nextclade_data/data/datasets/$d/versions/$NEW/files/tree.json
  aws s3 cp s3://nextstrain-staging/nextclade_sars-cov-2${EXTENSION[$d]}.json s3://nextstrain-data/nextclade_sars-cov-2${EXTENSION[$d]}.json
done
```

Edit CHANGELOG.md

## Regenerate escape.json after data update

Baseline needs to match the dataset reference, for XBB we don't ignore reversions as BA.2.86 has real reversions compared to XBB that aren't artefacts.

```bash
python scripts/generate_escape_calc_json.py \
  --baseline "BA.2" \
  --ignore-reversions \
  --wuhan defaults/wuhan_spike_translation.fasta \
  --output profiles/clades/21L/pathogen_json/escape.json
```

## Regenerate ace2.json after data update

```bash
python scripts/generate_ace2_json.py \
  --wuhan-spike defaults/wuhan_spike_translation.fasta \
  --reference-spike profiles/clades/XBB/spike_translation.fasta  \
  --output profiles/clades/XBB/pathogen_json/ace2.json
```

## Updating virus properties

```sh
#!/bin/bash
cd ~/code
for d in BA.2.86 wuhan-hu-1/{orfs,proteins}; do
  F=./nextclade_data/data/nextstrain/sars-cov-2/$d/pathogen.json
  jq --slurpfile v ./nextclade_data_workflows/sars-cov-2/virus_properties.json '.mutLabels.nucMutLabelMap=$v[0].nucMutLabelMap' $F | sponge $F
done
```
