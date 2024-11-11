import os
import pandas as pd
import multiprocessing as mp


## Files
wgs_file = 'processed_data/whole_genome_screened_sample_ids.tsv'
sample_file = 'raw_data/Cosmic_Sample_v100_GRCh38.tsv'
mutations_file = 'processed_data/Cosmic_GenomeScreensMutant_v100_GRCh38_filtered_snv_tss.tsv'
output_file = mutations_file[:-4] + '_wgs.tsv']


# Helper function to find samples that have been whole genome screened
def get_whole_genome_ids(chunk):
    return set(chunk[chunk['WHOLE_GENOME_SCREEN'] == 'y']['COSMIC_SAMPLE_ID'].unique())

# Get the COSS (unique COSMIC sample IDs) of the samples that were whole genome screened
if os.path.exists(wgs_file): # check if file already exists
	whole_genome_df = pd.read_csv(coss_genome_screened_file, sep='\t')
	whole_genome_ids = set(whole_genome_df.iloc[:, 0])
else: # if the file doesn't already exist then create it
	whole_genome_ids = set()
	with mp.Pool(mp.cpu_count()) as pool:
		chunks = pd.read_csv(sample_file, sep='\t', chunksize=10000)
		results = pool.map(get_whole_genome_ids, chunks)
		for result in results:
			whole_genome_ids.update(result)
			whole_genome_df = pd.DataFrame(whole_genome_ids, columns=['COSMIC_SAMPLE_ID'])
			whole_genome_df.to_csv(wgs_file, sep='\t', index=False)
del whole_genome_df


# Helper function to check if mutation is from whole genome screened samples
def filter_for_whole_genome(chunk):
    return chunk[chunk['COSMIC_SAMPLE_ID'].isin(whole_genome_ids)]

# Filter for mutations from whole genome screening
with pd.option_context('mode.chained_assignment', None):  # Suppress pandas warning for chained assignments
	with mp.Pool(mp.cpu_count()) as pool:
		mutation_chunks = pd.read_csv(mutations_file, sep='\t', chunksize=5000)
		filtered_chunks = pool.map(filter_for_whole_genome, mutation_chunks)
		if filtered_chunks:
			df = pd.concat(filtered_chunks)
			df.to_csv(output_file, sep='\t', index=False, header=True)

