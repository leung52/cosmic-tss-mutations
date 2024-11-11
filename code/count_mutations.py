import pandas as pd


## Files
tss_file = 'processed_data/hg38_tss_ranges_reduced.tsv'
snv_file = 'processed_data/Cosmic_GenomeScreensMutant_v100_GRCh38_filtered_snv_tss.tsv'
output_file = 'processed_data/wgs_mutation_count.tsv'


# Prepare gene coordinate ranges as dictionary
tss_df = pd.read_csv(tss_file, sep='\t')

tss_dict = {}
for chrom, chrom_group in tss_df.groupby('Chromosome'):
	tss_dict[chrom] = {}
	for strand, group in chrom_group.groupby('Strand'):
		tss_ranges = group[['TSS_Frame_Start', 'TSS_Frame_End']].values.tolist()
		sorted_ranges = sorted(tss_ranges, key=lambda x: x[0])
		tss_dict[chrom][strand] = sorted_ranges

# Finding relative TSS range
if tss_df['Strand'].iloc[0] == '+':
	upstream = tss_df['Transcriptional_Start_Site'].iloc[0] - tss_df['TSS_Frame_Start'].iloc[0]
	downstream = tss_df['TSS_Frame_End'].iloc[0] - tss_df['Transcriptional_Start_Site'].iloc[0]
else:
	upstream = tss_df['TSS_Frame_End'].iloc[0] - tss_df['Transcriptional_Start_Site'].iloc[0]
	downstream =  tss_df['Transcriptional_Start_Site'].iloc[0] - tss_df['TSS_Frame_Start'].iloc[0]

# Setting up mutation counts for each position in TSS window
count_dict = {}
for i in range(-upstream, downstream + 1):
	count_dict[i] = {
		'A>G': 0, 'A>C': 0, 'A>T': 0,
		'G>A': 0, 'G>C': 0, 'G>T': 0,
		'C>A': 0, 'C>G': 0, 'C>T': 0,
		'T>A': 0, 'T>G': 0, 'T>C': 0
	}

# Prevent repeat counting (key: COSS; value: COSV)
unique_dict = {}


# Nucleotide complement dict
comp_nuc = {
	'A': 'T',
	'T': 'A',
	'C': 'G',
	'G': 'C'
}


# Filter COSMIC tsv for SNV within TSS windows

for snv_df in pd.read_csv(snv_file, sep='\t', chunksize=5000):
	for _, row in snv_df.iterrows():
		coss = row['COSMIC_SAMPLE_ID']
		cosv = row['GENOMIC_MUTATION_ID']
		if coss not in unique_dict:
			unique_dict[coss] = [cosv]
		else:
			if cosv in unique_dict[coss]:
				continue
			unique_dict[coss].append(cosv)

		chrom = 'chr' + str(row['CHROMOSOME'])
		pos = row['GENOME_START']
		ref = row['GENOMIC_WT_ALLELE']
		alt = row['GENOMIC_MUT_ALLELE']

		if chrom in tss_dict:
			for (start, end) in tss_dict[chrom]['+']:
				if start <= pos <= end:
					rel_pos = pos - start - upstream
					mutation_type = ref + '>' + alt
					count_dict[rel_pos][mutation_type] += 1
			for (start, end) in tss_dict[chrom]['-']:
				if start <= pos <= end:
					rel_pos = (pos - start - downstream) * -1
					mutation_type = comp_nuc[ref] + '>' + comp_nuc[alt]
					count_dict[rel_pos][mutation_type] += 1


count_df = pd.DataFrame.from_dict(count_dict, orient='index').reset_index()
count_df.rename(columns={'index': 'Position'}, inplace=True)
count_df.to_csv(output_file, sep='\t', index=False)
