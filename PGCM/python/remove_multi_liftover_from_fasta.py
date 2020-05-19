"""
Filter a fasta file produced by MAKER
liftover (est2genome). In case the same
record id appears more than once, only
keep the one with best AED. If equal - 
take longest. If still equal - choose
one at random.
"""

from __future__ import print_function
from Bio import SeqIO
import sys

in_fasta = sys.argv[1]
out_fasta = sys.argv[2]

genes = {}
# go over fasta and select records by name
records = []
with open(in_fasta) as f:
  for rec in SeqIO.parse(f, "fasta"):
    gene_name = rec.id
    gene_len = len(rec.seq)
    aed = float(rec.description.split()[2].split(':')[1])
    if gene_name not in genes or aed < genes[gene_name][0] or (aed == genes[gene_name][0] and gene_len > genes[gene_name][1]):
      genes[gene_name] = (aed, gene_len, rec)

with open(out_fasta, 'w') as fo:
  SeqIO.write([r[2] for r in genes.values()], fo, "fasta")
