# Nextclade reference tree workflow for monkeypox

## Usage

```bash
snakemake -j all -pr -F
```

### Visualize results

View results with:

```bash
nextstrain view auspice/
```

## Maintenance

### Updating for new clades

- [ ] Update each `config/{build}/clades.tsv` with new clades
- [ ] Add new clades to color ordering
- [ ] Check that clades look good, exclude problematic sequences as necessary

### Creating a new dataset version

```sh
#!/bin/bash
OLD="2022-11-03T12:00:00Z"
NEW="2023-01-26T12:00:00Z"

B1="hMPXV_B1/references/pseudo_ON563414"
MPXV="MPXV/references/ancestral"
HPXV="hMPXV/references/NC_063383.1"

declare -A EXTENSION=( ["$B1"]="b1" ["$MPXV"]="MPXV" ["$HPXV"]="hpxv1" ) 

for d in $B1 $MPXV $HPXV; do
  rm -rf ~/code/nextclade_data/data/datasets/$d/versions/$NEW
  cp -pr ~/code/nextclade_data/data/datasets/$d/versions/$OLD ~/code/nextclade_data/data/datasets/$d/versions/$NEW
  sed -i "s/$OLD/$NEW/g" ~/code/nextclade_data/data/datasets/$d/versions/$NEW/files/tag.json
  aws s3 cp s3://nextstrain-staging/nextclade_monkeypox_${EXTENSION[$d]}.json - | gzcat >~/code/nextclade_data/data/datasets/$d/versions/$NEW/files/tree.json 
done
```

Edit CHANGELOG.md

## Configuration

Builds differ in paths, relevant configs are pulled in through lookup.

Configuration takes place in `config/config.yml` by default.
The analysis pipeline is contained in `workflow/snakemake_rule/core.smk`.
This can be read top-to-bottom, each rule specifies its file inputs and output and pulls its parameters from `config`.
There is little redirection and each rule should be able to be reasoned with on its own.

### Data use

We gratefully acknowledge the authors, originating and submitting laboratories of the genetic
sequences and metadata for sharing their work. Please note that although data generators have
generously shared data in an open fashion, that does not mean there should be free license to
publish on this data. Data generators should be cited where possible and collaborations should be
sought in some circumstances. Please try to avoid scooping someone else's work. Reach out if
uncertain.

## Installation

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.
Please choose the installation method for your operating system which uses Docker, as currently a pre-release version of Nextalign is required which we've baked into the `--image` argument to `nextstrain build` above.
