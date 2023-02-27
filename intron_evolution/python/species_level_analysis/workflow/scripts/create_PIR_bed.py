"""
Creates a BED file with random
Proximal Intergenic Regions (PIRs).
Takes a GFF3 file with genes and
introns, then:
1. Generate distribution of intron lengths
2. Create flanking regions with length X around genes
3. Create PIRs, randomly chosen from the flanking
   regions, with lengths corresponding to the
   intron lengths
4. Write to BED
"""

import sys
import random

genes_flank_bed = sys.argv[1]
introns_bed = sys.argv[2]
flank_size = int(sys.argv[3])
out_pir_bed = sys.argv[4]

# Parse GFF - create gene flanks and intron lengths
intron_len = []
flanks = []

with open(genes_flank_bed) as f:
  for line in f:
    fields = line.strip().split('\t')
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    flanks.append((chrom,start,end))

with open(introns_bed) as f:
  for line in f:
    fields = line.strip().split('\t')
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    l = end - start
    l = min(l, flank_size)
    intron_len.append(l)

pirs = []
for il in intron_len:
  pir = None
  while not pir:
    rand_flank = random.choice(flanks)
    rand_flank_chrom = rand_flank[0]
    rand_flank_start = rand_flank[1]
    rand_flank_end = rand_flank[2]
    rand_flank_len = rand_flank_end - rand_flank_start
    # if flank is too short - pick another
    if rand_flank_len < il:
      continue
    # choose pir start position in flank
    pir_start = random.randint(rand_flank_start, rand_flank_end-il)
    pir_end = pir_start + il
    pir = (rand_flank_chrom, str(pir_start), str(pir_end))
    pirs.append(pir)

with open(out_pir_bed, 'w') as fo:
  for pir in pirs:
    print('\t'.join(pir), file=fo)
