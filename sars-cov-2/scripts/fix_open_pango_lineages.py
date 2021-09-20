import argparse
from collections import defaultdict
import numpy as np
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="generate priorities files based on hash of strain name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", type=str, required=True)
    parser.add_argument("--designations", type = str, required=True)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t')
    designations = pd.read_csv(args.designations, sep=',',names=['strain','designation'])

    metadata['pango_designated'] = metadata.join(designations.set_index('strain'), on='strain')['designation']
    metadata['pango_lineage'] = metadata['pango_designated'].fillna(metadata['pango_lineage'])
    metadata.to_csv(args.output, sep='\t', index=False)
    

