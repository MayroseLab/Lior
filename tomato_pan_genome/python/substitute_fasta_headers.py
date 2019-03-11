"""
Perform string replacements in fasta headers based on a conversion table
"""

import sys

in_fasta = sys.argv[1]
in_tsv = sys.argv[2]
out_fasta = sys.argv[3]

convert_dict = {}
with open(in_tsv) as f:
  for line in f:
    from_str, to_str = line.strip().split('\t')
    convert_dict[from_str] = to_str

with open(in_fasta) as f, open(out_fasta,'w') as fo:
  for line in f:
    line = line.strip()
    if line.startswith('>'):
      for s in convert_dict:
        line = line.replace(s, convert_dict[s])
    print(line, file=fo)
