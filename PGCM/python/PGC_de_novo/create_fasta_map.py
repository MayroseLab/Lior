"""
Create a map file for use as input to
map_fasta_ids, based on a matching gff
"""

import sys

in_gff = sys.argv[1]

with open(in_gff) as f:
  for line in f:
    if line.startswith('#'):
      continue
    fields = line.strip().split('\t')
    if fields[2] == "mRNA":
      attr = {k.split('=')[0]: k.split('=')[1] for k in fields[8].split(';')[:-1]}
      if attr['ID'] == attr['Name']:
        print("%s\t%s" %(attr['Alias'], attr['ID']))
      else:
        print("%s\t%s" %(attr['Name'], attr['ID']))
