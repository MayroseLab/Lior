"""
Takes a gff file and for genes with
multiple mRNAs, only keeps the longest
(largest sum of exon lengths). Prints
out a new gff.
"""

from __future__ import print_function
import gffutils
import sys

in_gff = sys.argv[1]
out_gff = sys.argv[2]
gene_mrna_out = out_gff + '.gene_to_mRNA'

db_path = in_gff + '.db'
gff_db = gffutils.create_db(in_gff, db_path, force=True, merge_strategy="create_unique")
gff = gffutils.FeatureDB(db_path)

with open(out_gff, 'w') as fo, open(gene_mrna_out,'w') as fo2:
  for feature in gff.all_features():
    if feature.featuretype not in {'gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'}:
      print(str(feature), file=fo)
      continue
    if feature.featuretype != 'gene':	# mRNA and exon features
      continue
    gene = feature
    print(str(gene), file=fo)
    mrnas = list(gff.children(gene, featuretype='mRNA'))
    if len(mrnas) == 0:
      continue
    longest_transcript = mrnas[0]
    max_len = 0
    for mrna in mrnas:
      exons = list(gff.children(mrna, featuretype='exon'))
      total_exons_len = sum([exon.end - exon.start for exon in exons])
      if total_exons_len > max_len:
        max_len = total_exons_len
        longest_transcript = mrna
    print(str(longest_transcript), file=fo)
    print("%s\t%s" %(gene['ID'][0], longest_transcript['ID'][0]),file=fo2)
    for x in gff.children(longest_transcript):
      print(str(x), file=fo)