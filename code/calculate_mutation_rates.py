import pandas as pd


## Files:
nucleotide_counts_file = 'processed_data/nucleotide_count.tsv'
mutation_counts_file = 'processed_data/ncv_mutation_count.tsv'
output_file = 'processed_data/ncv_mutation_rates.tsv'

tss_nucleotide_counts = pd.read_csv(nucleotide_counts_file, sep='\t')
mutation_counts = pd.read_csv(mutation_counts_file, sep='\t')

mutation_types = ['A>G', 'A>C', 'A>T', 'G>A', 'G>C', 'G>T', 'C>A', 'C>G', 'C>T', 'T>A', 'T>G', 'T>C']

for mutation in mutation_types:
        mutation_counts[mutation] = mutation_counts[mutation] / tss_nucleotide_counts[mutation[0]]
        mutation_counts[mutation] = mutation_counts[mutation].replace([float('inf'), -float('inf')], pd.NA)

mutation_counts.to_csv(output_file, sep='\t', index=False)
