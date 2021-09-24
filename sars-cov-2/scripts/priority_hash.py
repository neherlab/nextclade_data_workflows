import argparse
from collections import defaultdict
import numpy as np
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="generate priorities files based on hash of strain name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", type=str, required=True, help="strain name txt")
    parser.add_argument("--seed", type = float, default=0, required=False, help="seed for stability")
    parser.add_argument("--output", type=str, required=True, help="tsv file with the priorities")
    args = parser.parse_args()

    priorities = pd.read_csv(args.input, sep='\t',header=None)

    priorities['priority'] = priorities.apply(lambda x: hash(x[0] + str(args.seed)) % 1000000, axis=1)

    priorities.to_csv(args.output, sep='\t', header=False, index=False)
    

