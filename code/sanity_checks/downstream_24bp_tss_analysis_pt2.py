import pandas as pd


snv_df = pd.read_csv('processed_data/downstream_24bp_tss_cosmic_snv.tsv', sep='\t')

counts = {
	'G>A': {},
	'G>T': {},
	'G>C': {}
}

for _, row in snv_df.iterrows():
	chrom = row['CHROMOSOME']
	pos = row['GENOME_START']
	mut = row['MUTATION_TYPE']
	coord = chrom + '.' + str(pos)
	if coord not in counts[mut]:
		counts[mut][coord] = 1
	else:
		counts[mut][coord] += 1

counts_df = pd.DataFrame(
	[(mut, coord, count) for mut, coords in counts.items() for coord, count in coords.items()],
	columns=['MUTATION_TYPE', 'COORDINATE', 'COUNT']
)

counts_df = counts_df.sort_values(by='COUNT', ascending=False)

counts_df.to_csv('processed_data/downstream_24bp_tss_mutation_counts.tsv', sep='\t', index=False)
