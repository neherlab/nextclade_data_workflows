## Copying produced files from scicore into dataset

```bash
scp -rC roemer0001@login-transfer.scicore.unibas.ch:~/nextclade_gisaid/sars-cov-2/auspice/nextclade/auspice.json sars-cov-2/references/MN908947/versions/2022-07-26T12:00:00Z/files/tree.json
```

## Deplying from staging to production

Replace date as appropriate:

```bash
cd auspice/nextclade
rm nextclade_*
cp auspice.json nextclade_sars-cov-2.json
cp nextclade_sars-cov-2.json nextclade_sars-cov-2_2022-07-26.json
cp auspice_raw_root-sequence.json nextclade_sars-cov-2_2022-07-26_root-sequence.json
nextstrain deploy s3://nextstrain-data nextclade_sars-cov-2*.json
```
