# Nextclade reference tree workflow for SARS-CoV-2

## Running the workflow

```sh
snakemake -c10 --profile profiles/clades -F
```

## Creating a new dataset version

```sh
#!/bin/bash
OLD="2023-04-18T12:00:00Z"
NEW="2023-05-10T12:00:00Z"

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
