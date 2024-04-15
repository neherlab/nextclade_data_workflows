# 

"""declare -A EXTENSION=( ["wuhan-hu-1/orfs"]="" ["wuhan-hu-1/proteins"]="_proteins" ["XBB"]="_XBB" ["BA.2.86"]="_BA.2.86" ["BA.2"]="_BA.2" )

for d in $SC2 $SC221L; do
    aws s3 cp s3://nextstrain-staging/nextclade_sars-cov-2${EXTENSION[$d]}.json - | gzcat >~/code/nextclade_data/data/datasets/$d/versions/$NEW/files/tree.json
    aws s3 cp s3://nextstrain-staging/nextclade_sars-cov-2${EXTENSION[$d]}.json s3://nextstrain-data/nextclade_sars-cov-2${EXTENSION[$d]}.json
done
"""

# New dataset involves:
# 1. Copying the tree.json file to the ~/code/nextclade_data/data/datasets/sars-cov-2/$PATH/tree.json
# 2. PATH = [wuhan-hu-1/orfs, wuhan-hu-1/proteins, XBB, BA.2.86, BA.2]
# 3. Mapping where the tree json is at various aws paths:

import os

path_to_url = {
    "wuhan-hu-1/orfs": "",
    "wuhan-hu-1/proteins": "_proteins",
    "XBB": "_XBB",
    "BA.2.86": "_BA.2.86",
    "BA.2": "_BA.2"
}

def full_path(path_part, filename="tree.json"):
    return f"~/code/nextclade_data/data/nextstrain/sars-cov-2/{path_part}/{filename}"

def full_from_url(url_part):
    return f"s3://nextstrain-staging/nextclade_sars-cov-2{url_part}.json"

def full_to_url(url_part):
    return f"s3://nextstrain-data/nextclade_sars-cov-2{url_part}.json"

for path_part, url_part in path_to_url.items():
    # Copy to nextclade_data
    os.system(f"aws s3 cp {full_from_url(url_part)} - | gzcat > {full_path(path_part)}")
    # Update README.md, prepend "defaults/README.md" at the path
    changelog_path = full_path(path_part, "CHANGELOG.md")
    os.system(f"cp defaults/CHANGELOG.md {changelog_path}_temp")
    os.system(f"echo '' >> {changelog_path}_temp")
    os.system(f"cat {changelog_path} >> {changelog_path}_temp")
    os.system(f"mv {changelog_path}_temp {changelog_path}")
    # Publish new tree
    # os.system(f"aws s3 cp {full_from_url(url_part)} {full_to_url(url_part)}")

