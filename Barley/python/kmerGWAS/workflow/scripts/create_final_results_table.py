import sys
import pandas as pd

pval_tsv = sys.argv[1]
direct_mapping_sam = sys.argv[2]
projection_tsvs = sys.argv[3:-1]
out_tsv = sys.argv[-1]

df_list = []

pval_df = pd.read_csv(pval_tsv, sep='\t', index_col=0)
df_list.append(pval_df)

direct_mapping_df = pd.read_csv(direct_mapping_sam, sep='\t', comment='@', usecols=[0,2,3], index_col=0, names=['k-mer_ID', 'Chr','Pos'])
if direct_mapping_df.empty:
  direct_mapping_df = pd.DataFrame(columns=['k-mer_ID', 'Chr','Pos'])
df_list.append(direct_mapping_df)

for proj_tsv in projection_tsvs:
  proj_df = pd.read_csv(proj_tsv, sep='\t', index_col=0).iloc[:,[2,3]]
  df_list.append(proj_df)

final_df = pd.concat(df_list, axis=1)
final_df.to_csv(out_tsv, sep='\t')
