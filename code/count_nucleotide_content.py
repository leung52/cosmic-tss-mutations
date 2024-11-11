import pandas as pd
from collections import defaultdict


## Files
tss_file = 'processed_data/hg38_tss_ranges_reduced.tsv'
fasta_file = 'processed_data/tss_sequences.fa'
output_file = 'processed_data/nucleotide_count.tsv'


# Find relative TSS range
tss_df = pd.read_csv(tss_file, sep='\t', nrows=2)

if tss_df['Strand'].iloc[0] == '+':
	upstream = tss_df['Transcriptional_Start_Site'].iloc[0] - tss_df['TSS_Frame_Start'].iloc[0]
	downstream = tss_df['TSS_Frame_End'].iloc[0] - tss_df['Transcriptional_Start_Site'].iloc[0]
else:
	upstream = tss_df['TSS_Frame_End'].iloc[0] - tss_df['Transcriptional_Start_Site'].iloc[0]
	downstream =  tss_df['Transcriptional_Start_Site'].iloc[0] - tss_df['TSS_Frame_Start'].iloc[0]


# Get nucleotide sequences
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


# Count nucleotides at each position in all sequences
count_dict = defaultdict(lambda: defaultdict(int))
for sequence in sequences:
        if len(sequence) == upstream + downstream + 1:
                for i, nucleotide in enumerate(sequence):
                        position = i - upstream
                        count_dict[position][nucleotide] += 1

count_dict = {pos: dict(counts) for pos, counts in count_dict.items()}

df = pd.DataFrame.from_dict(count_dict, orient='index')

df.index.name = 'Position'
df.reset_index(inplace=True)

df.to_csv(output_file, sep='\t', index=False)
