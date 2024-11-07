import os
import pandas as pd
import multiprocessing as mp


def get_whole_genome_ids(chunk):
    return set(chunk[chunk['WHOLE_GENOME_SCREEN'] == 'y']['COSMIC_SAMPLE_ID'].unique())

def filter_for_whole_genome(chunk):
    return chunk[chunk['COSMIC_SAMPLE_ID'].isin(whole_genome_ids)]

# Get a set of COSMIC sample IDs that were whole genome screened
id_file = 'processed_data/whole_genome_screened_sample_ids.tsv'

if os.path.exists(id_file):
	whole_genome_df = pd.read_csv(id_file, sep='\t')
	whole_genome_ids = set(whole_genome_df.iloc[:, 0])
else:
	sample_file = 'raw_data/Cosmic_Sample_v100_GRCh38.tsv'
	whole_genome_ids = set()

	with mp.Pool(mp.cpu_count()) as pool:
		chunks = pd.read_csv(sample_file, sep='\t', chunksize=10000)
		results = pool.map(get_whole_genome_ids, chunks)
		for result in results:
			whole_genome_ids.update(result)

			whole_genome_df = pd.DataFrame(whole_genome_ids, columns=['COSMIC_SAMPLE_ID'])
			whole_genome_df.to_csv(id_file, sep='\t', index=False)
del whole_genome_df

# Filter for mutations from whole genome screening
mutations_file = 'processed_data/Cosmic_NonCodingVariants_v100_GRCh38_filtered_SNV_and_TSS.tsv'

with pd.option_context('mode.chained_assignment', None):  # Suppress pandas warning for chained assignments
	with mp.Pool(mp.cpu_count()) as pool:
		mutation_chunks = pd.read_csv(mutations_file, sep='\t', chunksize=5000)
		filtered_chunks = pool.map(filter_for_whole_genome, mutation_chunks)

		if filtered_chunks:
			df = pd.concat(filtered_chunks)
			df.to_csv('processed_data/cosmic_whole_genome_screen_mutations.tsv', sep='\t', index=False, header=True)

