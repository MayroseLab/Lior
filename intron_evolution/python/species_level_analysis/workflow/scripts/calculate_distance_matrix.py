"""
Calculate pairwise distances
between species, based on
intron length vectors.
We use the KS statistic as
the distance metric.
"""

import sys
import os
from itertools import combinations
import pandas as pd
from scipy import stats

def gff_to_intron_lengths(gff):
  headers = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
  gff_df = pd.read_csv(gff, sep='\t', names=headers, comment='#')
  gff_df = gff_df.query('type == "intron"')
  return gff_df['end'] - gff_df['start'] + 1

if __name__ == "__main__":

  out_matrix = sys.argv[1]
  gff_list = sys.argv[2:]

  print('Parsing GFF files...')
  intron_lengths = {}
  for gff in gff_list:
    sp = os.path.basename(os.path.dirname(gff))
    print(sp)
    intron_lengths[sp] = gff_to_intron_lengths(gff)

  print('Calculating pairwise distances...')
  ks_similarity = []
  for pair in combinations(intron_lengths.keys(),2):
    sp1, sp2 = pair
    print(sp1, sp2)
    ks = stats.kstest(intron_lengths[sp1], intron_lengths[sp2]).statistic
    res = pd.Series([sp1,sp2,ks])
    ks_similarity.append(res)
    res_rev = pd.Series([sp2,sp1,ks])
    ks_similarity.append(res_rev)
  
  for sp in intron_lengths.keys():
    dist_to_self = pd.Series([sp,sp,0])
    ks_similarity.append(dist_to_self)

  print('Generating distance matrix...')
  ks_similarity_df = pd.concat(ks_similarity, axis=1).T
  ks_similarity_df.columns = ['species1','species2','KS']
  ks_similarity_df.set_index('species1', inplace=True)
  matrix = pd.pivot_table(ks_similarity_df, values='KS', index='species1', columns='species2')
  matrix.to_csv(out_matrix, sep='\t')
