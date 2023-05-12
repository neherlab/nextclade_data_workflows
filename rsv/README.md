# RSV dataset workflow

## Running the workflow

```bash
snakemake --profile profiles
```

## Testing datasets

Serve the dataset directory:

```bash
npx serve --cors  rsv/output/rsv_a/references/EPI_ISL_412866/versions/2023-02-03T12:00:00Z/files
```

Now open Nextclade in the browser and use the following URL:

```url
https://clades.nextstrain.org/?dataset-url=http://localhost:3000
```
