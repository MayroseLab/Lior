"""
Filters a SAM file of contigs
mapping to a reference genome:
* Remove unmapped contigs
* Remove mappings with MAPQ < x
* Remove duplicate mappings:
  - Choose the one with best MAPQ
  - If more than one exists, choose
    the one with most M in CIGAR
"""

import argparse
import re
from io import StringIO
import pandas as pd

def get_sam_header(sam_p):
  header = ''
  with open(sam_p) as f:
    line = f.readline()
    while line.startswith('@'):
      header += line
      line = f.readline()
  return header.strip()

def sam_to_df(sam_p):
  sam_headers = ["QNAM", "FLAG", "RNAM", "POS", "MAPQ", "CIGAR", "RNEX", "PNEX", "TLEN", "SEQ", "QUAL"]
  df = pd.read_csv(sam_p, sep='\t', comment='@', usecols=range(11), names=sam_headers)
  tags_df = pd.read_csv(sam_p, sep='|', comment='@', names=['ALL'])
  tags_df['TAGS'] = tags_df['ALL'].apply(lambda s: '\t'.join(s.split('\t')[11:]))
  tags_df.drop('ALL', axis=1, inplace=True)
  return pd.concat([df, tags_df], axis=1)

def list_to_pairs(l):
  it = iter(l)
  return list(zip(it,it))

def count_cigar_matches(cigar_str):
  cigar_pairs = list_to_pairs(re.split(r'([SHMID])', cigar_str))
  return sum([int(p[0]) for p in cigar_pairs if p[1] == 'M'])

def choose_mapping(sub_df):
  sub_df['CIGAR_m'] = sub_df['CIGAR'].apply(count_cigar_matches)
  sort_sub_df = sub_df.sort_values(by=['MAPQ','CIGAR_m'], ascending=False)
  sort_sub_df.drop('CIGAR_m', axis=1, inplace=True)
  return sort_sub_df.iloc[0]

def filter_sam(sam_df, min_mapq=0):
  # remove low quality mappings
  filt_sam_df = sam_df.query('MAPQ >= @min_mapq')
  # if no significant contigs remain
  if filt_sam_df.empty:
    return filt_sam_df
  # otherwise, remove duplicates
  filt_sam_df = filt_sam_df.groupby("QNAM").apply(choose_mapping)

  return filt_sam_df

def write_sam(sam_header, sam_df, out_sam):
  with open(out_sam, 'w') as fo:
    print(sam_header, file=fo)
    sio = StringIO()
    sam_df.to_csv(sio, sep='\t', mode='a', index=False, header=False)
    sio.seek(0)
    sam_str = sio.read().replace('"','')
    print(sam_str, file=fo)

if __name__ == "__main__":
  
  parser = argparse.ArgumentParser(description='Filter contig mapping SAM')
  parser.add_argument('in_sam', help='input SAM')
  parser.add_argument('min_mapq', type=int, help='Minimum MAPQ to keep')
  parser.add_argument('out_sam', help='output SAM')
  args = parser.parse_args()
  
  sam_header = get_sam_header(args.in_sam)
  sam_df = sam_to_df(args.in_sam)
  filt_sam_df = filter_sam(sam_df, args.min_mapq)
  write_sam(sam_header, filt_sam_df, args.out_sam)
