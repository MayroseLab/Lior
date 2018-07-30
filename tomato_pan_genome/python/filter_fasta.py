import argparse
from Bio import SeqIO
import re
import sys

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('in_fasta')
  parser.add_argument('-out_fasta')
  parser.add_argument('-e', nargs = '*')
  parser.add_argument('-v', nargs = '*')
  parser.add_argument('-lt', type=int)
  parser.add_argument('-gt', type=int)
  args = parser.parse_args()
  if args.gt and args.lt:
    assert args.gt < args.lt, "lt option must be larger than gt option"

  if args.e:
    positive_regex_list = [re.compile(r) for r in args.e]
  else:
    positive_regex_list = []
  if args.v:
    negative_regex_list = [re.compile(r) for r in args.v]
  else:
    negative_regex_list = []
  out_recs = []
  for rec in SeqIO.parse(args.in_fasta,'fasta'):
    if not all(r.search(rec.description) for r in positive_regex_list):
      continue
    if any(r.search(rec.description) for r in negative_regex_list):
      continue
    if args.gt and  not len(rec.seq) > args.gt:
      continue
    if args.lt and  not len(rec.seq) < args.lt:
      continue
    out_recs.append(rec)
    
  if not args.out_fasta or args.out_fasta == '-':
    out_handle = sys.stdout
  else:
    out_handle = open(args.out_fasta,'w')
  SeqIO.write(out_recs,out_handle,'fasta')
  out_handle.close()
