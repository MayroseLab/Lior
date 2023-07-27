"""
Given a GFF3 file containing intron features,
calculate various stats regarding the introns
and print to a stats file
"""

import sys
import os
import numpy as np
import pandas as pd
import scipy


### FUNCTIONS

def attributes_to_dict(attr):
  return {a.split('=')[0]: a.split('=')[1] for a in attr.split(';')}

def attributes_to_column(row, attr):
  attr_dict = attributes_to_dict(row['attributes'])
  if attr in attr_dict:
    return attr_dict[attr]
  else:
    return None
  
def gff_to_df(gff, type=None):
  headers = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
  gff_df = pd.read_csv(gff, sep='\t', names=headers, comment='#')
  if type:
    gff_df = gff_df.query('type == @type')
  gff_df['ID'] = gff_df.apply(attributes_to_column, axis=1, args=('ID',))
  gff_df['Parent'] = gff_df.apply(attributes_to_column, axis=1, args=('Parent',))
  return gff_df
  
def kde(v, steps=1000):
  kde = scipy.stats.gaussian_kde(v)
  vmin = v.min()
  vmax = v.max()
  step = (vmax-vmin)/steps
  x = np.arange(vmin,vmax,step)
  y = kde.evaluate(x)
  return x,y

def kde_peaks(x, y, frac=0.1, min_dist=20):
  """
  frac: height fraction of peak from previous highest peak to be keps
  min_dist: minimal distance on X axis between peaks
  """
  # find all peaks
  yrange = y.max() - y.min()
  peaks = scipy.signal.find_peaks(y, distance=min_dist, prominence=yrange/20)
  peaks_x = x[peaks[0]]
  peaks_y = y[peaks[0]]
  # convert peaks to (x,y pairs)
  peaks_pairs = list(zip(peaks_x, peaks_y))
  # sort peaks by descending y
  peaks_pairs_sort = sorted(peaks_pairs, key=lambda p: p[1], reverse=True)
  # iterate on peaks until peak y is smaller than frac of the previous peak
  i = 1
  while i < len(peaks_pairs_sort):
    if peaks_pairs_sort[i][1]/peaks_pairs_sort[i-1][1] >= frac:
      i += 1
    else:
      break
  return peaks_pairs_sort[:i]

def calc_stats(gff_df, column, name, query=None):
  
  stats_list = ['Min', 'Max','Mean', 'STD', 'Q10','Q25','Q50','Q75','Q90','Modes', 'M1_x', 'M1_y', 'M2_x', 'M2_y', 'Introns_count', 'Total_intron_length', 'Transcripts_count', 'Total_transcript_length', 'Transcripts_containing_introns', 'Mean_per_transcript', 'Total_intron_fraction', 'Mean_intron_fraction', 'Dataset']

  # per-transcript stats
  trans_df = gff_df.query('type == "mRNA"')
  df = gff_df.query('type == "intron"')

  df.sort_values(by=['seqid','start'], inplace=True)
  df_plus = df.query('strand == "+"')
  df_plus['intron_index'] = df_plus.groupby('attributes').cumcount()
  df_minus = df.query('strand == "-"')
  df_minus['intron_index'] = df_minus.groupby('attributes').cumcount(ascending=False)
  df = pd.concat([df_plus, df_minus])
  df['intron_index'] = df['intron_index'] + 1

  if query:
    df = df.query(query)
  if df.shape[0] == 0:
    stats = pd.Series([np.nan]*len(stats_list), index=stats_list)
    stats = pd.DataFrame(stats).transpose()
    return stats

  trans_count = trans_df.shape[0]
  trans_with_introns = df['Parent'].nunique()
  introns_per_trans = df.groupby('Parent')['intron_index'].max()
  single_exon_trans = set(trans_df['ID']) - set(introns_per_trans.index)
  single_exon_trans_introns = pd.Series(0, index=single_exon_trans)
  introns_per_trans = pd.concat([introns_per_trans, single_exon_trans_introns])
  mean_introns_per_trans = introns_per_trans.mean()

  trans_len = trans_df['length']
  trans_len.index = trans_df['ID']
  trans_len.name = 'transcript_length' 
  total_trans_len = trans_len.sum()
  
  intron_len_per_trans = df.groupby('Parent')[column].sum()
  single_exon_trans_introns_len = pd.Series(0, index=single_exon_trans)
  intron_len_per_trans = pd.concat([intron_len_per_trans, single_exon_trans_introns_len])
  intron_len_per_trans.name = column
  intron_len_per_trans = pd.concat([intron_len_per_trans, trans_len], axis=1)
  intron_len_per_trans['intron_fraction'] = intron_len_per_trans[column]/intron_len_per_trans['transcript_length']
  mean_intron_frac = intron_len_per_trans['intron_fraction'].mean()

  # introns stats

  vec = df[column]
  q10 = vec.quantile(0.1)
  q25 = vec.quantile(0.25)
  q50 = vec.quantile(0.5)
  q75 = vec.quantile(0.75)
  q90 = vec.quantile(0.9)
  min_ = vec.min()
  max_ = vec.max()
  mean = vec.mean()
  std = vec.std()
  
  kde_x, kde_y = kde(vec)
  peaks = kde_peaks(kde_x, kde_y)

  n_modes = len(peaks)
  if n_modes > 0:
    m1_x, m1_y = peaks[0]
  else:
    m1_x, m1_y = (np.nan, np.nan)
  if n_modes > 1:
    m2_x, m2_y = peaks[1]
  else:
    m2_x, m2_y = (np.nan, np.nan)

  count = df.shape[0]
  mean_per_transcript = df.groupby('Parent')['intron_index'].max().mean()

  total_len = vec.sum()
  total_intron_fraction = total_len/total_trans_len

  stats = pd.Series([min_, max_, mean, std, q10, q25, q50, q75, q90, n_modes, m1_x, m1_y, m2_x, m2_y, count, total_len,  trans_count, total_trans_len, trans_with_introns, mean_introns_per_trans, total_intron_fraction, mean_intron_frac, name], index=stats_list)

  stats = pd.DataFrame(stats).transpose()
  return stats

### MAIN

if __name__ == "__main__":

  in_gff = sys.argv[1]
  out_stats = sys.argv[2]
  species = sys.argv[3]

  gff_df = gff_to_df(in_gff)
  gff_df['length'] = gff_df['end'] - gff_df['start'] + 1
  gff_df['log_length'] = np.log10(gff_df['length'])

  all_introns_stats = calc_stats(gff_df, 'length', 'all')
  all_introns_log_stats = calc_stats(gff_df, 'log_length', 'all_log')
  first_introns_stats = calc_stats(gff_df, 'length', 'first', query='intron_index == 1')
  first_introns_log_stats = calc_stats(gff_df, 'log_length', 'first_log', query='intron_index == 1')
  nonfirst_introns_stats = calc_stats(gff_df, 'length', 'nonfirst', query='intron_index > 1')
  nonfirst_introns_log_stats = calc_stats(gff_df, 'log_length', 'nonfirst_log', query='intron_index > 1')

  stats_df = pd.concat([all_introns_stats,all_introns_log_stats,first_introns_stats,first_introns_log_stats,nonfirst_introns_stats,nonfirst_introns_log_stats])

  stats_df.index = [species]*stats_df.shape[0]

  stats_df.to_csv(out_stats, sep='\t') 
