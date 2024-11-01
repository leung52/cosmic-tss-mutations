# cosmic-tss-mutations

## Run Files order
Filter and count COSMIC mutations
```
python code/filter_cosmic_ncv_data.py
python code/count_mutations.py
```
Counts nucleotides
```
python code/make_bed.py
bedtools getfasta -fi raw_data/Homo_sapiens.GRCh38.dna.toplevel.fa -bed processed_data/tss_range_coordinates.bed -fo processed_data/tss_sequences.fa -s
python code/count_nucleotide_content.py
```
Calculate and graph mutation rate
```
python code/calculate_mutation_rates.py
python code/graph.py
```


## Steps taken:
1. Generate genome coordinate ranges for each TSS (-10 to 40)
        - run: 'tss_windows.py'
        - output: 'hg38_TSS_window.txt'
2. Filter non coding variants for mutants within -10 < TSS < 40
        - run: 'filter_cosmic_vcf.py'
        - output: 'Filtered_Cosmic_NonCodingVariants_v100_GRCh38.vcf'
3. Count number of mutations of each type for each position relative to TSS
        - run: 'count_mutations.py'
        - output: 'NonCodingVariants_Count.tsv'
4. Count number of each nucleotide at each position relative to TSS (exclude coding regions)
        - run: 'convert_to_BED.py'
        - output: 'tss_sequences.bed'
          - run: 'clean_up_hg38.py'
          - output: 'GCF_000001405.26_GRCh38_genomic_reduced.fna'
          - in bash: $ bedtools getfasta -fi GCF_000001405.26_GRCh38_genomic_reduced.fna -bed tss_sequences.bed -fo tss_sequences.fa
        - output: 'tss_sequences.fa'
          - run: 'count_nucleotides.py'
          - output: 'tss_nucleotide_counts.tsv'
6. Calculate mutation rates
          - run: 'calculate_mutation_rate.py'
          - output: 'noncoding_mutation_rates.tsv'
7. Graph mutation rates
          - run: 'graph_mutations.py'
          -output: 'mutation_rates_interactive.html'
