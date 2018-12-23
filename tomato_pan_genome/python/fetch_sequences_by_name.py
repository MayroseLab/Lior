from __future__ import print_function
from Bio import SeqIO
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_fasta',help='Input fasta')
parser.add_argument('-o','--output_fasta', help='output fasta')
parser.add_argument('-l','--list',help='file containing list of headers')
parser.add_argument('-v', default=False, action='store_true', help='Select all that do *not* match')
parser.add_argument('-p', '--partial_match', default=False, action='store_true', help='Allow partial matches')
args = parser.parse_args()

with open(args.list) as f:
  headers = set([l.strip() for l in f])

with open(args.input_fasta) as f, open(args.output_fasta,'w') as fo:
  for record in SeqIO.parse(f, "fasta"):
    match = False
    rec_header = record.name
    if args.partial_match:
      for h in headers:
        if h in rec_header:
          match = True
          break
    else:
      if rec_header in headers:
        match = True
    if (match and not args.v) or (not match and args.v):
      print(">%s" % rec_header, file=fo)
      print(record.seq, file=fo)

