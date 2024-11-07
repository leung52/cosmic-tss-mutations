import pandas as pd



df = pd.read_csv('processed_data/downstream_24bp_tss_cosmic_snv.tsv', sep='\t')


df = df.query('GENOME_START == 76736877')
print(df)



df.to_csv('processed_data/cosmic_snv_chr17_76736877.tsv', sep='\t', index=False)

