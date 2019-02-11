"""
Treat proteins sequences including a stop signal '*'
by breaking the sequence on '*' and taking longest
part, thereby ensuring no *'s are present.
"""

from __future__ import print_function
from Bio import SeqIO
import sys

in_fasta = sys.argv[1]
out_fasta = sys.argv[2]

with open(in_fasta, 'r') as f, open(out_fasta,'w') as fo:
  for record in SeqIO.parse(f, "fasta"):
    print(">%s" % record.description, file=fo)
    if '*' in record.seq:
      seq_breaks = record.seq.split('*')
      longest = max(seq_breaks, key=len)
      print(longest, file=fo)
    else:
      print(record.seq, file=fo)
