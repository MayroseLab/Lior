"""
Takes a gff file and for genes with
multiple mRNAs, only keeps the longest
(largest sum of exon lengths). Prints
out a new gff.
"""

from __future__ import print_function
import gffutils
import sys
from os import remove

in_gff = sys.argv[1]
out_gff = sys.argv[2]

db_path = "tmp.sqlite3"
gff_db = gffutils.create_db(in_gff, db_path, force=True, merge_strategy="create_unique")
gff = gffutils.FeatureDB(db_path)

with open(out_gff, 'w') as fo:
  for feature in gff.all_features():
    if feature.featuretype not in {'gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'}:
      print(str(feature), file=fo)
      continue
    if feature.featuretype != 'gene':	# mRNA and exon features
      continue
    gene = feature
    print(str(gene), file=fo)
    mrnas = list(gff.children(gene, featuretype='mRNA'))
    longest_transcript = mrnas[0]
    max_len = 0
    for mrna in mrnas:
      exons = list(gff.children(mrna, featuretype='exon'))
      total_exons_len = sum([exon.end - exon.start for exon in exons])
      if total_exons_len > max_len:
        max_len = total_exons_len
        longest_transcript = mrna
    print(str(longest_transcript), file=fo)
    for x in gff.children(longest_transcript):
      print(str(x), file=fo)
        
remove(db_path)
