"""
Combine results from multiple SGSGeneLoss
.excov files into one PAV matrix:
	sample_1|sample_2|...|sample_n
gene_1	1       |1       |...|0
gene_2	0       |1       |...|1
...
gene_m	1       |0       |...|1

The script takes a list of .excov files
(one per sample) sepearated by spaces.
File names should be:
<sample name>.all.PAV
The last argument is the output tsv.

"""

import pandas as pd
import sys
import os

in_files = sys.argv[1:-2]
ref_name = sys.argv[-2]
out_tsv = sys.argv[-1]

samples_pav = []
for f in in_files:
  sample_name = os.path.basename(f).split('.')[0]
  df = pd.read_csv(f)
  pav_s = df['is_lost'].map({'PRESENT': 1, 'LOST': 0})
  pav_s.name = sample_name
  pav_s.index = df['ID']
  pav_s.index.rename('gene', inplace=True)
  samples_pav.append(pav_s)

pav_df = pd.concat(samples_pav, axis=1)
pav_df[ref_name] = pd.Series(pav_df.index, index=pav_df.index).str.startswith('PanGene').map({False: 1, True: 0})
pav_df.to_csv(out_tsv, sep='\t')
