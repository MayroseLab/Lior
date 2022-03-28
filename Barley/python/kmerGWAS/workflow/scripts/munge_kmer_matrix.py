"""
Read k-mer counts matrix from stream,
filter, convert to presence-absence,
and extract stats.
The expected input from SDTIN are as
created by KMC matrixer:
kmer    count1    count2    ...    count_n
No header line, all tab-separated.
"""

import sys
import argparse


if __name__ == "__main__"

  parser = argparse.ArgumentParser(description='Read k-mer counts matrix from stream')
  
  parser.add_argument('--min_maf', '-m', type=float, default=0.05, help='Min minor allele frequency (0-1) to keep k-mers')
  parser.add_argument('--acc_list', '-a', default=None, help='File containing accession names list (one per line)')
  parser.add_argument('--acc_n', '-n', type=int, default=None, help='Number of accessions in matrix')
  parser.add_argument('out_CNV', help='Output path for k-mer CNV matrix')
  parser.add_argument('out_PAV', help='Output path for k-mer PAV matrix')
  parser.add_argument('out_stats', help='Output path k-mer stats file')

  args = parser.parse_args()

  # Accession names

  if args.acc_list:
    with open(args.acc_list) as f:
      acc_names = [l.strip() for l in f.readlines()]
  elif args.acc_n:
    acc_names = ["Acc_%s" % i for i in range(args.acc_n)]
  else:
    sys.exit('Must specify accession names list OR number of accessions')
  acc_num = len(acc_names)
  fields_valid_n = acc_num + 1

  # Data structures for storing stats

  kmer_occup = {c: 0 for c in range(1, acc_num+1)}

  acc_kmer_occup = {acc: {c: 0 for c in range(1, acc_num+1)} for acc in acc_names}

  pa_patterns = {}

  # Read data from STDIN

  for line in sys.stdin:
    fields = line.strip().split()
    kmer = fields.pop(0)
    kmer_pa = counts_to_pa(fields)
 
