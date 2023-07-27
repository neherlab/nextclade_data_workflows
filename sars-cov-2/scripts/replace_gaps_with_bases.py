#%%
from Bio import SeqIO

def replace_gaps_with_bases(seq1, seq2):
    '''This function replaces gaps in seq1 with corresponding bases from seq2.'''

    # Check if both sequences are of same length
    if len(seq1) != len(seq2):
        print("Error: Both sequences should be of the same length")
        return

    new_seq1 = []
    for base1, base2 in zip(seq1, seq2):
        if base1 == '-':
            new_seq1.append(base2)
        else:
            new_seq1.append(base1)

    return ''.join(new_seq1)

#%%

# Reading two sequences from fasta file
seq_records = list(SeqIO.parse("xbb.fasta", "fasta"))
seq1 = str(seq_records[0].seq)
seq2 = str(seq_records[1].seq)

# Replace gaps in seq1 with corresponding bases from seq2
new_seq2 = replace_gaps_with_bases(seq2, seq1)

#%%
new_seq2_list = list(new_seq2)
new_seq2_list[15939-1] = 'C'  # Python indexing is 0-based, so we subtract 1
new_seq2 = ''.join(new_seq2_list)

# Write new sequence to fasta file
with open("xbb_in_wuhan_coordinates.fasta", "w") as f:
    f.write(f">{seq_records[1].id}\n{new_seq2}")

# %%
