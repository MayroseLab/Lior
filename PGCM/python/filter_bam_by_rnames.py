"""
filters a bam file based on a list
of RNAME values, so only reads with
RNAME or and RNEXT in the list are
kept. Unlike samtools, the header is
adjusted - all other SQ entries are
removed. Automatically indexes the
output bam.
"""

import sys
import pysam

in_bam = sys.argv[1]
in_list = sys.argv[2]
out_bam = sys.argv[3]

bam = pysam.AlignmentFile(in_bam)
required_rnames = set([l.strip() for l in open(in_list)])

print("Adjusting bam header")
bam_header = bam.header.to_dict()
required_rids = [ bam.get_tid(rname) for rname in required_rnames ]
new_header = bam_header.copy()
new_header["SQ"] = [new_header["SQ"][i] for i in required_rids]
new_bam = pysam.AlignmentFile(out_bam, 'wb', header=new_header)

print("Starting to write new bam")
c = 0
for read in bam.fetch():
  if bam.get_reference_name(read.reference_id) not in required_rnames or bam.get_reference_name(read.next_reference_id)  not in required_rnames:
    continue
  old_rid = read.reference_id
  old_rname = bam.get_reference_name(old_rid)
  new_rid = new_bam.get_tid(old_rname)
  read.reference_id = new_rid
  old_rid_next = read.next_reference_id
  old_rname_next = bam.get_reference_name(old_rid_next)
  new_rid_next = new_bam.get_tid(old_rname_next)
  read.next_reference_id = new_rid_next
  new_bam.write(read)
  c += 1
  if c % 1000000 == 0:
    print("%s reads written..." % c)

bam.close()
new_bam.close()
pysam.index(out_bam)
