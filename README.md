## Checking new tree

1. Download generated files into nextclade data workflow repo:

    ```bash
    scp -rC roemer0001@login-transfer.scicore.unibas.ch:~/nextclade_data_workflows/sars-cov-2/output output
    ```

1. Plug them into nextclade.org advanced view.
1. Filter to new nodes and check that:
    - clades are clean
    - no big outliers
1. Check `tag.json` is up to date (ideally update in `profiles/tag.json` for posterity)
1. Check `qc.json` does not regress (ideally update in `profiles/qc.json` for posterity) [beware, codons are 0 indexed]
1. Potentially run `scripts/common_stops.py` and `scripts/common_frameshifts.py` to add new stops/frameshifts that have become more common to `qc.json`

## Committing to data repo

1. Go to nextclade_data_workflow repo
1. Checkout branch, open PR to master
1. Copy output from workflow repo to data repo

    ```bash
    cp -r output/sars-cov-2/references/MN908947/versions/  ../../nextclade_data/data/datasets/sars-cov-2/references/MN908947/versions
    ```

1. Update `changelog.md`
1. Get Ivan to review
1. Merge into master

## Release process

Follow release guidelines as outlined here: https://github.com/nextstrain/nextclade_data#dataset-release-process
