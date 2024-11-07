import pandas as pd

# File with the genomic TSS information (must have 'Transcription_Start_Site' as a column name)
input_file = 'raw_data/Fantom5_CAGE_hg38_gencode.v40_TSS_EIB_Updated -200.txt'

output_file = 'processed_data/hg38_tss_ranges.tsv'

# Helper function to calculate the genomic coordinates for the TSS range
def calculate_frames(row):

	upstream_read = 10
	downstream_read = 40

	# Positive strand
	if row['Strand'] == '+':
        	frame_start = row['Transcriptional_Start_Site'] - upstream_read
        	frame_end = row['Transcriptional_Start_Site'] + downstream_read
	# Negative strand
	else:
        	frame_start = row['Transcriptional_Start_Site'] - downstream_read
        	frame_end = row['Transcriptional_Start_Site'] + upstream_read

	return pd.Series({'TSS_Frame_Start': frame_start, 'TSS_Frame_End': frame_end})


tss_coordinates_df = pd.read_csv(input_file, sep=r'\s+')
tss_coordinates_df[['TSS_Frame_Start', 'TSS_Frame_End']] = tss_coordinates_df.apply(calculate_frames, axis=1)
tss_reading_frame_df = tss_coordinates_df[['Chromosome', 'Transcriptional_Start_Site', 'TSS_Frame_Start', 'TSS_Frame_End', 'Strand']].copy()
tss_reading_frame_df.to_csv(output_file, sep='\t', index=False)
