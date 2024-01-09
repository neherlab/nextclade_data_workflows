import argparse
import json
from typing import Any, Dict

def main(args):
    # Your main logic here
    # Example: print(args.labeled_muts, args.reversions, args.escape, args.ace2, args.output)
    # Load each file
    
    result: dict[str, Any] = json.load(open(args.base_json, "r"))

    labeled = json.load(open(args.labeled_muts, "r"))

    from Bio import SeqIO
    wuhan = str(SeqIO.read(args.wuhan, "fasta").seq)
    reference = str(SeqIO.read(args.reference, "fasta").seq)
    
    for i in range(len(wuhan)):
        if wuhan[i] != reference[i]:
            labeled["nucMutLabelMap"][f"{i+1}{wuhan[i]}"] = ["rev"]
    
    # Add labeled mutations
    result["mutLabels"] = { **labeled }

    escape = json.load(open(args.escape, "r"))
    ace2 = json.load(open(args.ace2, "r"))
    
    if escape != {} or ace2 != {}:
        result["phenotypeData"] = []
    
    if escape != {}:
        result["phenotypeData"].append(escape)
    
    if ace2 != {}:
        result["phenotypeData"].append(ace2)

    

    extras = json.load(open(args.extras, "r"))
    result.update(extras)

    with open(args.output, "w") as f_out:
        json.dump(result, f_out, indent=2)
        f_out.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Assemble pathogen JSON.')
    parser.add_argument('--base-json', required=True, help='Path to base pathogen JSON file')
    parser.add_argument('--labeled-muts', required=True, help='Path to labeled mutations file')
    parser.add_argument('--wuhan', required=True, help='Path to wuhan fasta')
    parser.add_argument('--reference', required=True, help='Path to reference fasta')
    parser.add_argument('--escape', required=True, help='Path to escape file')
    parser.add_argument('--ace2', required=True, help='Path to ACE2 file')
    parser.add_argument('--extras', required=True, help='Path to extras file')
    parser.add_argument('--output', required=True, help='Path for output pathogen JSON file')

    args = parser.parse_args()
    main(args)