import pandas as pd


## Files
tss_file = 'processed_data/hg38_tss_ranges_reduced.tsv'
output_file = "processed_data/tss_range_coordinates.bed"

df = pd.read_csv(tss_file, sep='\t')

# Create BED formatted dataframe
bed_df = bed_df = pd.DataFrame({
        'chrom': df['Chromosome'].str.replace('chr', '', regex=False),
        'start': df['TSS_Frame_Start'] - 1,  # BED is 0-based
        'end': df['TSS_Frame_End'],
        'name': df['Chromosome'] + ':' + df['Transcriptional_Start_Site'].astype(str),
        'extra': '.',
        'strand': df['Strand']
})

bed_df.to_csv(output_file, sep='\t', header=False, index=False)
