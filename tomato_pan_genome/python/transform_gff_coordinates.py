"""
Takes a gff file, where record names have the form
X_start_end and converts to record name X with start
and end columns adjusted.
This is useful when merging split gff files.
"""
from __future__ import print_function
import sys

in_gff = sys.argv[1]
out_gff = sys.argv[2]

with open(in_gff) as f, open(out_gff,'w') as fo:
  for line in f:
    line = line.strip()
    if line.startswith('#'):
      print(line,file=fo)
      continue
    fields = line.split('\t')
    if len(fields) != 9:
      print(line,file=fo)
      continue
    if fields[2] == "contig":
      continue
    rec_name_with_coords = fields[0]
    rec_name_split = fields[0].split('_')
    start_on_real_rec = int(rec_name_split[-2])
    real_rec_name = '_'.join(rec_name_split[:-2])
    fields[0] = real_rec_name
    fields[3] = str(int(fields[3]) + start_on_real_rec)
    fields[4] = str(int(fields[4]) + start_on_real_rec)
    print('\t'.join(fields), file=fo)
