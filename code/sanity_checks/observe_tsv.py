import pandas as pd

# Define the file path
file_path = 'raw_data/Cosmic_NonCodingVariants_v100_GRCh38.tsv'

# Initialize counters
total_rows = 0
unique_values = [set(), set(), set()]

# Define the column you want to find unique values for
columns = ['COSMIC_SAMPLE_ID', 'GENOMIC_MUTATION_ID', 'GENOMIC_MUT_ALLELE']

# Use pandas to read in chunks
chunk_size = 500  # Define a chunk size
for chunk in pd.read_csv(file_path, sep='\t', chunksize=chunk_size, usecols=columns):
	total_rows += len(chunk)  # Add rows in each chunk to the total
	for column, unique_set in zip(columns, unique_values):
		unique_set.update(chunk[column].dropna().unique())  # Update unique values set

# Print results
print("Total number of data points:", total_rows)
for column, unique_set in zip(columns, unique_values):
	if len(unique_set) > 20:
		print(f"Number of unique values in {column}: {len(unique_set)}")
	else:
		print(f"{column} unique values: {unique_set}")

