# Take a fasta file and create a new one
# only containing coordinates given in bed file.
# This is similar to bedtools getfasta, but more robust.

from __future__ import print_function
from Bio import SeqIO
import sys

in_fasta = sys.argv[1]
in_bed = sys.argv[2]
out_fasta = sys.argv[3]

bed_d = {}
with open(in_bed) as f:
  for line in f:
    fields = line.strip().split('\t')
    if len(fields) >= 3:
      chrom, start, end = fields[:3]
      start = int(start) - 1
      end = int(end) - 1
      name = "%s_%s_%s" %(chrom, start, end)
    else:
      exit("Line %s is not valid bed." % line.strip())
    if len(fields) == 4:
      name = fields[4]
    if chrom not in bed_d:
      bed_d[chrom] = {}
    bed_d[chrom][name] = (start, end)

with open(in_fasta) as f, open(out_fasta,'w') as fo:
  for rec in SeqIO.parse(f,'fasta'):
    rec_id =  rec.id
    if rec_id not in bed_d:
      continue
    rec_seq = rec.seq
    for x in bed_d[rec_id]:
      print('>%s' % x, file=fo)
      print(rec_seq[bed_d[rec_id][x][0]:bed_d[rec_id][x][1]+1], file=fo)
   
