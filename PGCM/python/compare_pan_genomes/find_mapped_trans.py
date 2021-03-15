from __future__ import print_function, division
import sys
import os

trans_fasta = sys.argv[1]
per_sample_dir = sys.argv[2]
out_dir = sys.argv[3]
out_list = sys.argv[4]

# generate transcripts dict
mapped = {}
with open(trans_fasta) as f:
  for line in f:
    if line.startswith('>'):
      name = line.strip()[1:]
      mapped[name] = [''] * 4

if not os.path.isdir(out_dir):
  os.mkdir(out_dir)

for sample in os.listdir(per_sample_dir):
  sample_dir = os.path.join(per_sample_dir, sample)
  assembly_dir = [d for d in os.listdir(sample_dir) if d.startswith('RG_assembly')]
  if not assembly_dir:
    continue
  genome_fasta = os.path.join(sample_dir, assembly_dir[0],'ragtag_output/ragtag.scaffolds.fasta')
  out_paf = os.path.join(out_dir, sample + '.paf')
  os.system('minimap2 -x splice:hq -uf %s %s > %s' %(genome_fasta, trans_fasta, out_paf))

  with open(out_paf) as f:
    for line in f:
      fields = line.strip().split('\t')
      if int(fields[9])/int(fields[1])  >= 0.95:
        mapped[fields[0]] = [sample, fields[5], fields[7], fields[8]]

with open(out_list,'w') as fo:
  for trans in mapped:
    print('\t'.join([trans] + mapped[trans]), file=fo)
