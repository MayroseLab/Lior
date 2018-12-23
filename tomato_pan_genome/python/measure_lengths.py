from Bio import SeqIO
import sys

in_fasta = sys.argv[1]
with open(in_fasta) as f:
  for rec in SeqIO.parse(f,'fasta'):
    print("%s\t%s" %(rec.id, len(rec.seq)))
