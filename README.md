# cosmic-tss-mutations
## Steps taken:
1. Generate genome coordinate ranges for each TSS (-10 to 40)
        * run: 'tss_windows.py'
        * output: 'hg38_TSS_window.txt'
2. Filter non coding variants for mutants within -10 < TSS < 40
        - run: 'filter_cosmic_vcf.py'
        - output: 'Filtered_Cosmic_NonCodingVariants_v100_GRCh38.vcf'
3. Count number of mutations of each type for each position relative to TSS
        - run: 'count_mutations.py'
        - output: 'NonCodingVariants_Count.tsv'
4. Count number of each nucleotide at each position relative to TSS (exclude coding regions)
        - run: 'convert_to_BED.py'
        - output: ''
        - in bash: $
        - output: ''
5. 
