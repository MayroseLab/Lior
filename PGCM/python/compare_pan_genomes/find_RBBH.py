"""
This script takes two reciprocal blast6
files and detects reciprocal best blast
hits (RBBHs). It is assumed that blast
was run with -max_target_seqs 1.
Outputs a blast6 tsv containing only RBBHs.
"""

from __future__ import print_function
import sys

in_blast6_tsv1 = sys.argv[1]
in_blast6_tsv2 = sys.argv[2]
out_RBBH_tsv = sys.argv[3]

# parse blast outputs and find RBBHs
def parse_blast6_file(p):
  res = {}
  with open(p) as f:
    for line in f:
      fields = line.strip().split('\t')
      res[fields[0]] = fields
  return res

fwd_res = parse_blast6_file(in_blast6_tsv1)
rev_res = parse_blast6_file(in_blast6_tsv2)

with open(out_RBBH_tsv, 'w') as fo:
  for prot in fwd_res:
    prot_best_hit = fwd_res[prot][1]
    if prot_best_hit in rev_res and rev_res[prot_best_hit][1] == prot:
      print('\t'.join(fwd_res[prot]), file=fo)
