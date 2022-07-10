"""
The purpose of this script is to project
or convert coordinates from one accession
to another. Specifically, we can project
k-mer mapping coordinates from contigs of
one assembly to reference coordinates.
To do that, we need two SAM files:
1. k-mers mapping to contigs
2. contigs mapping to reference
Given these two inputs, k-mer mappings
are projected onto the reference using the
following logic:
If the k-mer is mapped within a contig region
which is matched to the reference - use the
resulting reference position.
Otherwise (soft/hard clips and insertions),
use reference coordinates of the closest matched
region.
The output is a simple TSV table with k-mer
positions in reference coordinates.
"""

import pandas as pd
from intervaltree import Interval, IntervalTree
import sys
import re

def sam_to_df(sam_p, header_pref=''):
  """
  Read a SAM file as pandas DF.
  Ignore any SAM tags.
  header_pref allows to add a
  prefix to column names (e.g.
  genome name)
  """
  sam_headers = ["QNAM", "RNAM", "POS", "CIGAR"]
  df = pd.read_csv(sam_p, sep='\t', comment='@', usecols=[0,2,3,5], names=sam_headers, dtype={'POS': int})
  if header_pref:
    df.columns = ["%s_%s" %(header_pref, col) for col in df.columns]
  return df

def list_to_pairs(l):
  """
  Helper function to convert a
  list into a list of pair tuples
  with consequent elements:
  [1,2,3,4] -> [(1,2),(3,4)]
  """
  it = iter(l)
  return list(zip(it,it))

def cigar_to_ivtree(cigar_str):
  """
  Takes a SAM CIGAR string and
  returns an interval tree with
  coordinate ranges (starting from
  1) to their CIGAR operations:
  10M3I12M -> {[1:10]: 'M', [11:13]: 'I', [14:26]: 'M'}
  """
  cigar_pairs = list_to_pairs(re.split(r'([SHMID])', cigar_str))
  cigar_ivt = IntervalTree()
  s = 1
  for p in cigar_pairs:
    if p[1] not in 'SHMID':
      continue
    l = int(p[0])
    cigar_ivt[s:s+l] = p[1]
    s += l
  return cigar_ivt

def project_kmer_position(pos_on_contig, contig_cigar_ivt, start_from=1):
  """
  Project contig coordinate
  to reference coordinate.
  start_from is the starting
  coordinate of the contig on
  the reference.
  """
  c = 1
  for iv in sorted(contig_cigar_ivt):
    if iv.data == 'D':
      continue
    if iv.end < pos_on_contig + 1:
      if iv.data in 'SHM':
        c += (iv.end - iv.begin)
    else:
      if iv.data in 'SHM':
        res = start_from + pos_on_contig - iv.begin + c
      elif iv.data == 'I':
        res = start_from + pos_on_contig - iv.begin + c - (iv.end - iv.begin)/2
      return (int(res), iv.data)

if __name__ == "__main__":

  kmer_sam = sys.argv[1]
  contigs_sam = sys.argv[2]
  acc_name = sys.argv[3]
  out_tsv = sys.argv[4]

  # Read k-mers and contigs SAM
  kmer_sam_df = sam_to_df(kmer_sam, 'kmer')
  contigs_sam_df = sam_to_df(contigs_sam, 'contig')

  out_headers = ['k-mer_ID', 'contig', 'contig_Pos', 'Chr_%s' % acc_name, 'Pos_%s' % acc_name, 'Pos_type']
  if kmer_sam_df.empty or contigs_sam_df.empty:
    out_df = pd.DataFrame(columns=out_headers)
  else:
    # Join k-mers and contigs DFs based on contig name
    join_df = kmer_sam_df.merge(contigs_sam_df, how='inner', left_on='kmer_RNAM', right_on='contig_QNAM')
    # Convert CIGAR strings to interval trees
    join_df['contig_CIGAR_ivt'] = join_df['contig_CIGAR'].apply(cigar_to_ivtree)
    # Project k-mer coordinates
    join_df[['kmer_POS_project','kmer_POS_type']] = join_df.apply(lambda row: project_kmer_position(row['kmer_POS'], row['contig_CIGAR_ivt'], row['contig_POS']), axis=1, result_type="expand")
    # Create output DF: k-mer ID   contig   coordinate_on_contig   chromosome   coordinate_on_chromosome   CIGAR_operation
    # add acc_name to chromosome and coordinate_on_chromosome
    out_df = join_df[['kmer_QNAM', 'kmer_RNAM', 'kmer_POS', 'contig_RNAM', 'kmer_POS_project', 'kmer_POS_type']]
    out_df.columns = out_headers
  out_df.to_csv(out_tsv, sep='\t', index=False)
