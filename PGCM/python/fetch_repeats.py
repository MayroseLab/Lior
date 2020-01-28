"""
Scan a fasta file for stretches oflowercase
bases and output them to a new fasta as
separate records.
"""
from __future__ import print_function
from Bio import SeqIO
import sys
import re

in_fasta = sys.argv[1]

all_repeats = set()
for rec in SeqIO.parse(in_fasta,'fasta'):
  seq = rec.seq
  repeats = re.split(r'[ATGCN]+',str(seq))
  all_repeats = all_repeats.union(set(repeats))
rep_id = 1
for rep in all_repeats:
  if len(rep) < 2:
    continue
  print('>rep_%s' % rep_id)
  print(rep.upper())
  rep_id += 1
