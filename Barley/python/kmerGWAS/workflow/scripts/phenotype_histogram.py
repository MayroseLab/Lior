import sys
import pandas as pd
import plotly.express as px

in_tsv = sys.argv[1]
pheno_name = sys.argv[2]
out_html = sys.argv[3]

pheno_df = pd.read_csv(in_tsv, sep='\t')
n_accessions = pheno_df.shape[0]
plot_title = "%s (N=%s)" %(pheno_name.replace('_',' '), n_accessions)
max_val = pheno_df['phenotype_value'].max()
fig = px.histogram(pheno_df['phenotype_value'], title=plot_title)
fig.update_xaxes(range=[0, max_val])
fig.update_layout(showlegend=False)
fig.write_html(out_html, include_plotlyjs='cdn')

