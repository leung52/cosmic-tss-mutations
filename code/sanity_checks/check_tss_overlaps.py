import pandas as pd

# Load the TSS data
tss_df = pd.read_csv('processed_data/hg38_TSS_window.tsv', sep='\t')

# Create intervals for the TSS windows
tss_df['interval'] = tss_df.apply(lambda row: pd.Interval(left=row['TSS_Frame_Start'], right=row['TSS_Frame_End']), axis=1)

overlapping_rows = []

for chrom in tss_df['Chromosome'].unique():
	chrom_data = tss_df[tss_df['Chromosome'] == chrom]

	for i in range(len(chrom_data)):
		current_interval = chrom_data.iloc[i]['interval']
		overlaps = chrom_data[chrom_data['interval'].apply(lambda x: x.overlaps(current_interval)) & (chrom_data.index != chrom_data.index[i])]

		if not overlaps.empty:
			overlapping_rows.extend(overlaps.values.tolist())

# Create a DataFrame from the overlapping rows
overlapping_df = pd.DataFrame(overlapping_rows, columns=tss_df.columns)
overlapping_df = overlapping_df.drop(columns=['interval'])

overlapping_df.to_csv('processed_data/overlapping_tss_ranges.tsv', sep='\t', index=False)
