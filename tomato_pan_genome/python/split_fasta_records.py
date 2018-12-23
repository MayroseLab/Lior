# Take a fasta file and create a new one in which every record is
# split into chunks of a user-defined size. New record IDs are:
# ><original id>_<start>_<end>

from __future__ import print_function
from Bio import SeqIO
import sys

in_fasta = sys.argv[1]
out_fasta = sys.argv[2]
chunk_size = int(sys.argv[3])

with open(in_fasta) as f, open(out_fasta,'w') as fo:
  for rec in SeqIO.parse(f,'fasta'):
    rec_id =  rec.id
    rec_seq = rec.seq
    seq_len = len(rec_seq)
    if seq_len < chunk_size:
      chunk_id = "%s_0_%s" %(rec_id, seq_len)
      print(">%s" % chunk_id, file=fo)
      print(rec_seq, file=fo)
      continue
    for i in range(0,seq_len,chunk_size):
      start = i
      end = min(i + chunk_size, seq_len)
      chunk_id = "%s_%s_%s" %(rec_id,start,end)
      print(">%s" % chunk_id, file=fo)
      print(rec_seq[start:end], file=fo)
   
