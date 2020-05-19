"""
Takes a gff file produced by MAKER liftover
(est2genome) and in case an EST produced
multiple genes, only keeps one.
Prefers genes with lowest AED, if equal
takes longest, if still equal, choose one
at random
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

# go over all genes and create genes dict:
# {gene_name: {ID: (AED,length)}}
# select best gene
genes = {}
with open(out_gff, 'w') as fo:
  for gene in gff.features_of_type('gene'):
    gene_name = gene["Name"][0]
    gene_id = gene["ID"][0]
    gene_len = gene.end - gene.start
    mRNA = next(gff.children(gene, featuretype='mRNA'))
    aed = float(mRNA["_AED"][0])
    if gene_name not in genes or aed < genes[gene_name][1] or (aed == genes[gene_name][1] and gene_len > genes[gene_name][2]):
      genes[gene_name] = (gene_id, aed, gene_len)

# go over gff again and print only best genes
with open(out_gff, 'w') as fo:
  for feature in gff.all_features():
    if feature.featuretype not in {'gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'}:
      print(str(feature), file=fo)
      continue
    if feature.featuretype != 'gene':	# mRNA and exon features
      continue
    gene = feature
    if gene["ID"][0] == genes[gene["Name"][0]][0]:
      print(str(gene), file=fo)
      for c in gff.children(gene):
        print(str(c), file=fo)
        
remove(db_path)
