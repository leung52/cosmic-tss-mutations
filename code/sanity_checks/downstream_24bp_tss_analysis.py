import pandas as pd


## Prepare gene coordinate ranges as dictionary
tss_to_ignore = pd.read_csv('processed_data/deleted_tss_windows.tsv', sep='\t')
tss_df = pd.read_csv('processed_data/hg38_TSS_window.tsv', sep='\t').merge(tss_to_ignore, how='left', indicator=True)

tss_dict = {}

for chrom, chrom_group in tss_df.groupby('Chromosome'):
	tss_dict[chrom] = {}
	for strand, group in chrom_group.groupby('Strand'):
		tss_ranges = group[['TSS_Frame_Start', 'TSS_Frame_End']].values.tolist()
		sorted_ranges = sorted(tss_ranges, key=lambda x: x[0])
		tss_dict[chrom][strand] = sorted_ranges

# Nucleotide complement dict
comp_nuc = {
	'A': 'T',
	'T': 'A',
	'C': 'G',
	'G': 'C'
}


## Filter COSMIC tsv for SNV within TSS windows
snv_file = 'processed_data/Cosmic_NonCodingVariants_v100_GRCh38_filtered_SNV_and_TSS.tsv'
rows = []

for snv_df in pd.read_csv(snv_file, sep='\t', chunksize=5000):
	for _, row in snv_df.iterrows():
		chrom = 'chr' + str(row['CHROMOSOME'])
		pos = row['GENOME_START']
		ref = row['GENOMIC_WT_ALLELE']
		alt = row['GENOMIC_MUT_ALLELE']

		if chrom in tss_dict:
			for (start, end) in tss_dict[chrom]['+']:
				if start <= pos <= end:
					rel_pos = pos - start - 10
					mutation_type = ref + '>' + alt
					if rel_pos == 24 and ref == 'G':
						new_row = row.copy()
						new_row['MUTATION_TYPE'] = mutation_type
						new_row['STRAND'] = '+'
						rows.append(new_row)
			for (start, end) in tss_dict[chrom]['-']:
				if start <= pos <= end:
					rel_pos = (pos - start - 40) * -1
					mutation_type = comp_nuc[ref] + '>' + comp_nuc[alt]
					if rel_pos == 24 and comp_nuc[ref] == 'G':
						new_row = row.copy()
						new_row['MUTATION_TYPE'] = mutation_type
						new_row['STRAND'] = '-'
						rows.append(new_row)

result_df = pd.DataFrame(rows)
result_df.to_csv('processed_data/downstream_24bp_tss_cosmic_snv.tsv', sep='\t', index=False)

