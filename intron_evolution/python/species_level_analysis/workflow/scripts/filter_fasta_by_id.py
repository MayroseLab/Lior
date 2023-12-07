"""
Given a fasta file and a list of IDs,
only prints out records with rec.id
included in the list.
The list may be given as a file in the
second command line argument, or taken
from STDOIN
If the list contains a second field
(separate by tab), then it is used
for record name replacement.
"""

from Bio import SeqIO
import sys

# if ids file given
if len(sys.argv) == 3:
  inp = open(sys.argv[2]).readlines()
# if should be taken from STDIN
elif len(sys.argv) == 2:
  inp = sys.stdin
id_dict = {}
for l in inp:
  l = l.strip()
  fields = l.split('\t')
  if len(fields) == 0:
    continue
  elif len(fields) == 1:
    id_dict[fields[0]] = fields[0]
  else:
    id_dict[fields[0]] = fields[1]


rec_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))
for mrna_id in id_dict:
  if mrna_id in rec_dict:
    rec = rec_dict[mrna_id]
  elif mrna_id.split('.')[0] in rec_dict:
    rec = rec_dict[mrna_id.split('.')[0]]
  elif not mrna_id.endswith('.1') and mrna_id+'.1' in rec_dict:
    rec = rec_dict[mrna_id+'.1']
    rec.id = mrna_id
  else:
    continue
  rec.name = ''
  rec.description = ''
  print(rec.format('fasta').strip())
