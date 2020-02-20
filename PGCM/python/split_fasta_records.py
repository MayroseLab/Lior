"""
Split fasta records in one file into X files
"""

from __future__ import print_function, division
from Bio import SeqIO
import argparse
from itertools import zip_longest
from random import shuffle

parser = argparse.ArgumentParser()
parser.add_argument('in_fasta', help="input_fasta file")
parser.add_argument('out', help="Output fasta or output dir (in break mode)")
parser.add_argument('n_chunks', help="number of files to split into", type=int)
args = parser.parse_args()

def grouper(n, iterable, padvalue=None):
    "grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
    return zip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

n_records = len([1 for line in open(args.in_fasta) if line.startswith(">")])
records_per_chunk = max(n_records // args.n_chunks,1)
records = list(SeqIO.parse(args.in_fasta,'fasta'))
shuffle(records) # in case records are sorted by length
records_grouped = grouper(records_per_chunk,records)
i = 1
for g in records_grouped:
  g = [rec for rec in g if rec is not None] # avoid None records
  file_path = args.out + '/chunk%s.fasta' % i
  with open(file_path,'w') as fo:
    SeqIO.write(g, fo, "fasta")
  i += 1
