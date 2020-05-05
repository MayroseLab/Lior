"""
Given a SAM (or BAM) file, extract
sequences correspoding to alignment
inserts (I operations in CIGAR) larger
than the user-defined cutoff.
Print results to stdout in fasta format.
"""

import pysam
import sys

in_sam = sys.argv[1]
min_insert = int(sys.argv[2])

sam = pysam.AlignmentFile(in_sam)

for aln in sam:
  cigar_tuples = aln.cigartuples
  if not cigar_tuples:
    continue
  start = 0
  for ct in cigar_tuples:
    if ct[0] == 1 and ct[1] > min_insert and aln.query_sequence:	# insertion (I) operation
      insert_seq = aln.query_sequence[start:start+ct[1]]
      seq_name = "%s_%s_insert" %(aln.query_name, start)
      print('>%s' % seq_name)
      print(insert_seq)
    start += ct[1]
