"""
Search a fasta file for a given
regex, and output matches coordinates
in BED format
"""

import sys
import re
from Bio import SeqIO

in_fasta = sys.argv[1]
regex = sys.argv[2]
out_bed = sys.argv[3]


#regex = re.compile(regex.encode('unicode_escape'))
regex = re.compile(regex)
with open(out_bed, 'w') as fo:
  for rec in SeqIO.parse(in_fasta, 'fasta'):
    for m in regex.finditer(str(rec.seq)):
      l = [ rec.id, m.start(), m.start() + len(m.group()) ]
      l = [str(x) for x in l]
      print('\t'.join(l), file=fo)
