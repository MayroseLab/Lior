"""
This script iteratively maps chromosome-
or contig-level assemblies to a reference
genome, while adding novel sequences to
the reference at each iteration, thereby
creating a pan genome.
The input is a TSV file with header line:
sample	genome_fasta  annotation_gff (optional)
Genes on novel sequences will be extracted from
the gff3 files.
*** The order of genomes in the input TSV may change the results ***
"""

from __future__ import print_function
import os
import csv
import argparse

if __name__ == "__main__":

  # command line args
  parser = argparse.ArgumentParser()
  parser.add_argument('ref_fasta', help='Path to reference genome fasta')
  parser.add_argument('ref_gff', help='Path to reference annotation gff')
  parser.add_argument('in_tsv', help='Path to TSV file with samples data')
  parser.add_argument('out_dir', help='Path to output directory')
  parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use in Minimap')
  parser.add_argument('--min_len', default=1, type=int, help='Minimal sequence length to add to pan')
  args = parser.parse_args()

  # initialize
  pan_fasta = args.ref_fasta
  pan_gff = args.ref_gff

  # read TSV
  with open(args.in_tsv) as f:
    reader = csv.DictReader(f, delimiter = '\t')
    for row in reader:
      genome_name = row['sample']
      assembly_fasta = row['genome_fasta']
      genome_gff = None
      if 'annotation_gff' in row and row['annotation_gff']:
        genome_gff = row['annotation_gff']

      # run minimap2
      out_paf = os.path.join(args.out_dir, "%s_vs_pan.paf" % genome_name)
      os.system("minimap2 -x asm5 -L -t %s %s %s -o %s -c" % (cpus, pan_fasta, assembly_fasta, out_paf))

      # extract novel sequences
      script_dir = os.path.dirname(sys.argv[0])
      extract_script = os.path.join(script_dir, "extract_novel_regions_from_paf.py")
      novel_seq_fasta = os.path.join(args.out_dir, "%s_novel.fasta" % genome_name)
      if genome_gff:
        novel_gff = os.path.join(args.out_dir, "%s_novel.gff" % genome_name)
        os.system("python %s %s %s %s %s --in_gff %s --out_gff %s --genome_name %s" %(extract_script, out_paf, assembly_fasta, min_len, novel_seq_fasta, genome_gff, novel_gff, genome_name))
      else:
        os.system("python %s %s %s %s %s --genome_name %s" %(extract_script, out_paf, assembly_fasta, min_len, novel_seq_fasta, genome_name))

      # create new pan genome
      curr_pan = os.path.join(args.out_dir, "current_pan.fasta")
      os.system("cat %s %s > %s" %(pan_fasta, novel_seq_fasta, curr_pan + '.tmp'))
      os.replace(curr_pan + '.tmp', curr_pan)

      # if gff exists, add to pan genome annotation
      if genome_gff:
        curr_gff = os.path.join(args.out_dir, "current_pan.gff")
        os.system("cat %s %s > %s" %(pan_gff, novel_gff, curr_gff + '.tmp'))
        os.replace(curr_gff + '.tmp', curr_gff)

      # update pan
      pan_fasta = curr_pan
      if genome_gff:
        pan_gff = curr_gff

  # create final pan genome and gff
  final_pan = os.path.join(args.out_dir, "pan_genome.fasta")
  os.link(pan_fasta, final_pan)
  final_gff = os.path.join(args.out_dir, "pan_genes.gff")
  os.link(pan_gff, final_gff)
