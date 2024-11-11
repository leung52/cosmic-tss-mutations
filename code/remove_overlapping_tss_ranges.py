import pandas as pd



## Data files
tss_file = 'processed_data/hg38_tss_ranges.tsv'
output_file = 'processed_data/hg38_tss_ranges_reduced.tsv'

# Duplicate TSS for a gene
data = {
    "Chromosome": ["chr3", "chr3", "chr6", "chr8", "chr9", "chr11", "chr16", "chr19"],
    "Transcriptional_Start_Site": [184249695, 196338373, 37257653, 10839847, 128881931, 4697831, 15395754, 14835197],
    "TSS_Frame_Start": [184249685, 196338333, 37257613, 10839807, 128881891, 4697791, 15395744, 14835157],
    "TSS_Frame_End": [184249735, 196338383, 37257663, 10839857, 128881941, 4697841, 15395794, 14835207],
    "Strand": ["+", "-", "-", "-", "-", "-", "+", "-"]
}

tss_to_ignore = pd.DataFrame(data)

tss_df = pd.read_csv(tss_file, sep='\t')
tss_df = tss_df[~tss_df.set_index(['Chromosome', 'Transcriptional_Start_Site', 'TSS_Frame_Start', 'TSS_Frame_End', 'Strand'])
		.index.isin(tss_to_ignore.set_index(['Chromosome', 'Transcriptional_Start_Site', 'TSS_Frame_Start', 'TSS_Frame_End', 'Strand']).index)]

tss_df.to_csv(output_file, sep='\t', index=False)
