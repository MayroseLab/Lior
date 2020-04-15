"""
Split fasta records in one file into X files
Takes the number of chunks to split to and
calculates the required chunk size.
If needed (and allowed), will break records.
In such cases, the user may provide a bed file
that tells the script what regions should not be
broken.
"""

from __future__ import print_function, division
from Bio import SeqIO
import argparse
from intervaltree import IntervalTree

parser = argparse.ArgumentParser()
parser.add_argument('in_fasta', help="input_fasta file")
parser.add_argument('out', help="Output dir")
parser.add_argument('n_chunks', help="number of files to split into", type=int)
parser.add_argument('--no_breaks', action="store_true", default=False, help="Do not break records" )
parser.add_argument('--regions_bed', help="Bed file defining regions that shouldn't be broken")
args = parser.parse_args()

records = list(SeqIO.parse(args.in_fasta,'fasta'))
total_len = sum([len(rec) for rec in records])
len_per_chunk = total_len/args.n_chunks
out_chunks_bed = args.out + '/chunks.bed'

if args.regions_bed:
  regions = {}
  with open(args.regions_bed) as f:
    for line in f:
      rec_name, start, end = line.strip().split('\t')[:3]
      start = int(start)
      end = int(end)
      if rec_name not in regions:
        regions[rec_name] = IntervalTree()
      regions[rec_name][start:end] = True

curr_chunk_len = 0
curr_chunk_records = []
chunk_i = 1

for rec in records:
  rec_len = len(rec)
  if curr_chunk_len + rec_len <= len_per_chunk:
    curr_chunk_records.append(rec)
    curr_chunk_len += rec_len
  else:
    # first, print previous chunk (unless it's empty)
    if curr_chunk_records:
      file_path = args.out + '/chunk%s.fasta' % chunk_i
      SeqIO.write(curr_chunk_records, file_path, "fasta")
      # then reset chunk
      chunk_i += 1
      curr_chunk_len = 0
      curr_chunk_records = []
    # if the record is shorter than the max chunk size, add to new chunk
    if rec_len <= len_per_chunk:
      curr_chunk_len += rec_len
      curr_chunk_records.append(rec)
    # if the record is longer than the max chunk size, start breaking (if allowed)
    else:
      if args.no_breaks:
        curr_chunk_records.append(rec)
        curr_chunk_len = rec_len 
      else:
        k = 0
        while rec_len >= len_per_chunk:
          break_coord = k+int(len_per_chunk)
          if args.regions_bed:
            if rec.id in regions and regions[rec.id][break_coord]:
              break_coord = sorted(regions[rec.id][break_coord])[0].end + 1
          break_rec = rec[k:break_coord]
          break_rec.id += "__%s-%s" %(k,break_coord)
          file_path = args.out + '/chunk%s.fasta' % chunk_i
          SeqIO.write([break_rec], file_path, "fasta")
          chunk_i += 1
          rec_len -= break_coord - k
          k = break_coord
        if rec_len == 0:
          continue
        else:
          break_rec = rec[k:]
          break_rec.id += "__%s-%s" %(k,'end')
          file_path = args.out + '/chunk%s.fasta' % chunk_i
          SeqIO.write([break_rec], file_path, "fasta")
          chunk_i += 1
        curr_chunk_len = 0
        curr_chunk_records = []
file_path = args.out + '/chunk%s.fasta' % chunk_i
SeqIO.write(curr_chunk_records, file_path, "fasta")
