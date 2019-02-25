"""

"""

from __future__ import print_function, division
from Bio import SeqIO
import argparse
import os
from math import ceil

parser = argparse.ArgumentParser()
parser.add_argument('in_fasta', help="input_fasta file")
parser.add_argument('out', help="Output fasta or output dir (in break mode)")
parser.add_argument('--chunk_size', help="chunk size in bp", type=int)
parser.add_argument('--n_chunks', help="number of chunks to break_mode into", type=int)
parser.add_argument('--break_mode', default=False, action='store_true', help="write each chunk to a separate file")
args = parser.parse_args()

if not (args.chunk_size or args.n_chunks) or (args.chunk_size and args.n_chunks):
  exit("Must specify either chunk size or n chunks, but not both")
if args.break_mode and not os.path.isdir(args.out):
  exit("Output must be a directory when using break_mode mode")

# if number of chunks was given, get total sequence
# length and calculate required chunk size
if args.n_chunks:
  tot_len = 0
  with open(args.in_fasta) as f:
    for rec in SeqIO.parse(f,'fasta'):
      rec_seq = rec.seq
      tot_len += len(rec_seq)
  chunk_size = int(ceil(tot_len // args.n_chunks))
else:
  chunk_size = args.chunk_size

if not args.break_mode:	# single output file
  fo = open(args.out,'w')
with open(args.in_fasta) as f:
  for rec in SeqIO.parse(f,'fasta'):
    rec_id =  rec.id
    rec_seq = rec.seq
    seq_len = len(rec_seq)
    for i in range(0,seq_len,chunk_size):
      start = i
      end = min(i + chunk_size, seq_len)
      chunk_id = "%s_%s_%s" %(rec_id,start,end)
      chunk_seq = rec_seq[start:end]
      if args.break_mode:
        out_file = "%s/%s.fasta" %(args.out, chunk_id)
        fo = open(out_file,'w')
      print(">%s" % chunk_id, file=fo)
      print(chunk_seq, file=fo)
      if args.break_mode:
        fo.close()
if not args.break_mode:
  fo.close()
