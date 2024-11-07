import pandas as pd
from collections import defaultdict


## Files
fasta_file = 'processed_data/tss_sequences.fa'
output_file = 'processed_data/nucleotide_count.tsv'


sequences = []

with open(fasta_file, 'r') as f:
        sequence = ''
        for line in f:
                line = line.strip()
                if line.startswith('>'):
                        if sequence:
                                sequences.append(list(sequence.upper()))
                        sequence = ''
                else:
                        sequence += line
        if sequence:
                sequences.append(list(sequence.upper()))


# Initialize count dictionary for positions -10 to +40 using defaultdict
count_dict = defaultdict(lambda: defaultdict(int))

# Count nucleotides at each position in all sequences
for sequence in sequences:
        if len(sequence) == 51:
                for i, nucleotide in enumerate(sequence):
                        position = i - 10
                        count_dict[position][nucleotide] += 1

count_dict = {pos: dict(counts) for pos, counts in count_dict.items()}

df = pd.DataFrame.from_dict(count_dict, orient='index')

df.index.name = 'Position'
df.reset_index(inplace=True)

df.to_csv(output_file, sep='\t', index=False)
