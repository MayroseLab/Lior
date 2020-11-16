"""
Match two sets of proteins based on BLAST results,
using the maximum weight bipartite matching method.
Parameters:
1. set1 proteins fasta
2. set2 proteins fasta
3. set1 vs. set2 blast6 result
4. set2 vs. set1 blast6 result
5. Weight parameter - e.g. bitscore or pident
6. Min weight - minnimal weight value allowed in a match

Output:
TSV file with the matches pairs and their respective weight
"""

from __future__ import print_function, division
from networkx import bipartite, matching
from networkx.algorithms.bipartite.matching import minimum_weight_full_matching
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument('set1_fasta', help="Proteins set 1 fasta")
parser.add_argument('set2_fasta', help="Proteins set 2 fasta")
parser.add_argument('set1_vs_set2_blast', help="Blast6 result of set1 vs set2")
parser.add_argument('set2_vs_set1_blast', help="Blast6 result of set2 vs set1")
parser.add_argument('out_tsv', help="Path to output file")
parser.add_argument('--weight_param', default='bitscore', help="Blast column to use as edge weight" )
parser.add_argument('--min_weight', type=float, default=0, help="Minnimal weight value allowed in a match")
parser.add_argument('--set1_name', default="set1", help="Name of set1")
parser.add_argument('--set2_name', default="set2", help="Name of set2")
args = parser.parse_args()

# go over fasta files and assign an id to each protein
def assign_ids(fasta, start=0):
  i = start
  d = {}
  with open(fasta) as f:
    for line in f:
      if line.startswith('>'):
        d[i] = line.strip()[1:]
        i += 1
  return d

set1_proteins = assign_ids(args.set1_fasta)
set1_size = len(set1_proteins)
set2_proteins = assign_ids(args.set2_fasta, start=set1_size)
set2_size = len(set2_proteins)

# create bipartite graph
bg = bipartite.complete_bipartite_graph(set1_size,set2_size)

# parse blast6 files
def parse_blast6(blast6, weight_param):
  res = {}
  fieldnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
  with open(blast6) as f:
    reader = csv.DictReader(f, delimiter='\t', fieldnames=fieldnames)
    for row in reader:
      if (row['qseqid'], row['sseqid']) in res:
        w = float(res[(row['qseqid'], row['sseqid'])])
      else:
        w = 0
      res[(row['qseqid'], row['sseqid'])] = max(w,float(row[weight_param]))
  return res

print("Parsing BLAST files...")
set1_vs_set2_weights = parse_blast6(args.set1_vs_set2_blast, args.weight_param)
set2_vs_set1_weights = parse_blast6(args.set2_vs_set1_blast, args.weight_param)

print("Assigning edge weights...")
# go over all BG edges and assign weights
for edge in bg.edges:
  q = set1_proteins[edge[0]]
  s = set2_proteins[edge[1]]
  bidirectional_weight = 0
  if (q,s) in set1_vs_set2_weights:
    bidirectional_weight += set1_vs_set2_weights[(q,s)]
  if (s,q) in set2_vs_set1_weights:
    bidirectional_weight += set2_vs_set1_weights[(s,q)]
  bidirectional_weight /= 2
  bidirectional_weight *= -1 # b/c the function looks for the minimum weight matching
  bg.edges[edge]['weight'] = bidirectional_weight

print("Finding max weight matching...")
# find max weight matching
mw = minimum_weight_full_matching(bg)

print("Writing out results...")
# print out matches
with open(args.out_tsv, 'w') as fo:
  print("%s\t%s\t%s weight" %(args.set1_name, args.set2_name, args.weight_param), file=fo)
  for match in mw:
    source, target = match, mw[match]
    if target < source:
      # each match appears in both directions in the output dict
      continue
    match_weight = -1 * bg[source][target]['weight']
    if match_weight >= args.min_weight:
      q_name = set1_proteins[source]
      s_name = set2_proteins[target]
      print("%s\t%s\t%s" %(q_name, s_name, match_weight), file=fo)
