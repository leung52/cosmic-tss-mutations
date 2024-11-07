import pandas as pd
import multiprocessing as mp


## Data files
tss_file = 'processed_data/hg38_tss_ranges_reduced.tsv'
cosmic_ncv_file = 'raw_data/Cosmic_NonCodingVariants_v100_GRCh38.tsv'
output_file = 'processed_data/Cosmic_NonCodingVariants_v100_GRCh38_filtered_SNV_and_TSS.tsv'


# Prepare gene coordinate ranges as dictionary
tss_df = pd.read_csv(tss_file, sep='\t')
tss_dict = {}

for chrom, group in tss_df.groupby('Chromosome'):
	tss_ranges = group[['TSS_Frame_Start', 'TSS_Frame_End']].values.tolist()
	sorted_ranges = sorted(tss_ranges, key=lambda x: x[0])

	merged = []

	for start, end in sorted_ranges:
		# If merged is empty or there is no overlap, add the range
		if not merged or merged[-1][1] < start:
			merged.append([start, end])
		else:
			# There is an overlap, so we merge the current range with the last one
			merged[-1][1] = max(merged[-1][1], end)

	tss_dict[chrom] = merged


# Helper function for filtering
def filter_snv_in_tss_range(df):
	NUCLEOTIDES = {'A', 'T', 'C', 'G'}
	snv_df = df[df['GENOMIC_WT_ALLELE'].isin(NUCLEOTIDES) & df['GENOMIC_MUT_ALLELE'].isin(NUCLEOTIDES)]
	snv_in_tss_range_ls = []
	for _, row in snv_df.iterrows():
		chrom = f"chr{row['CHROMOSOME']}"
		pos = row['GENOME_START']
		# Check if chromosome exists in TSS dictionary and if position falls in any TSS interval
		if chrom in tss_dict:
			for start, end in tss_dict[chrom]:
				if start <= pos <= end:
					snv_in_tss_range_ls.append(row)
					break   # Stop checking intervals if a match is found
	if snv_in_tss_range_ls != []:
		return pd.DataFrame(snv_in_tss_range_ls)
	else:
		return pd.DataFrame(columns=df.columns)


# Filter COSMIC tsv for SNV within TSS windows
header_written = False

# Initialize multiprocessing pool
with mp.Pool(processes=mp.cpu_count()) as pool:
	with open(output_file, 'w') as f_out:
		# Process each chunk in parallel
		for filtered_df in pool.imap(filter_snv_in_tss_range, pd.read_csv(cosmic_file, sep='\t', chunksize=5000)):
			# Only write non-empty results
			if not filtered_df.empty:
				filtered_df.to_csv(f_out, sep='\t', index=False, header=not header_written)
				header_written = True  # Ensure header is only written once
