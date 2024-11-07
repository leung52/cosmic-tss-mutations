import plotly.graph_objects as go
import pandas as pd

## Files
data_file = 'processed_data/ncv_mutation_rates.tsv'
output_folder = 'graphs/noncoding_variants/'

## Define mutation type groups that separate graphs
mutation_groups = {
	"GC mutations": ['A>G', 'A>C', 'T>G', 'T>C'],
	"AT mutations": ['G>A', 'G>T', 'C>A', 'C>T'],
	"Transversions": ['C>G', 'G>C', 'T>A', 'A>T']
}



# Check if the folder exists, and create it if it doesn't
os.makedirs(output_folder, exist_ok=True)

df = pd.read_csv(data_file, sep='\t')
positions = df['Position']

# Loop through each mutation group and create a figure
for group_name, mutation_types in mutation_groups.items():
	fig = go.Figure()

	# Add a trace for each mutation type in the group
	for mutation_type in mutation_types:
		if mutation_type in df.columns:  # Only add if mutation type is in the data
			fig.add_trace(go.Scatter(
				x=positions,
				y=df[mutation_type],
				mode='lines',
				name=mutation_type
			))

	# Customize the layout
	fig.update_layout(
		title=f'Mutation Rates by Position ({group_name})',
		xaxis_title='Position',
		yaxis_title='Mutation Rate',
		legend_title='Mutation Types',
		template='simple_white',
		hovermode='x',
		xaxis=dict(showgrid=False),
		yaxis=dict(showgrid=False),
		font=dict(size=18)
	)

	fig.add_vline(x=0, line=dict(color="grey", width=2, dash="dash"))

	# Save each figure as an HTML file
	fig.write_html(output_folder + f"{group_name.replace(' ', '_').lower()}.html")
