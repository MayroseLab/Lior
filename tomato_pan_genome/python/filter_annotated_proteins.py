"""
This script is meant to be run on protein sets produced by an annotation process.
It takes in protein sequences and produces a filtered set based on some criteria
and cutoffs that can be set by the user.
"""

import argparse
import os
import sys
sys.path.append("../../queue_utilities/")
from queueUtils import send_commands_to_queue
import argparse
import logging
from pypette import Pipe, Job, BashJob


class FullPaths(argparse.Action):
  """Expand user- and relative-paths"""

  def __call__(self, parser, namespace, values, option_string=None):
    setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


class StreamToLogger(object):
  """
  Fake file-like stream object that redirects writes to a logger instance.
  """

  def __init__(self, logger, log_level=logging.INFO):
    self.logger = logger
    self.log_level = log_level
    self.linebuf = ''

  def write(self, buf):
    for line in buf.rstrip().splitlines():
      self.logger.log(self.log_level, line.rstrip())

  def flush(self):
    pass

def blastp():
  pass
def calc_c_score():
  pass
def calc_prot_coverage():
  pass
def calc_repeats_overlap():
  pass
def filter_proteins():
  pass

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('queue_conf', action=FullPaths)
  parser.add_argument('in_prot_fasta', help="Input protein sequences in fasta format", action=FullPaths)
  parser.add_argument('db', help="Pre-formated blast DB to run against", action=FullPaths)
  parser.add_argument('out_dir', action=FullPaths, help="Output directory")
  parser.add_argument('--in_prot_gff', default=None, action=FullPaths,help="A gff3 file containing genome coordinates for the proteins")
  parser.add_argument('--repeats_gff', default=None, action=FullPaths,help="A gff3 file containing repeat regions coordinates")
  parser.add_argument('--max_evalue', default=10.0, type=float, help="Ignore any blast hits with E-value larger than X")
  parser.add_argument('--c_score_cutoff', default=0.5, type=float, help="C score cutoff under which proteins are filtered")
  parser.add_argument('--prot_coverage_cutoff', default=0.5, type=float, help="Protein coverage cutoff under which proteins are filtered")
  parser.add_argument('--repeats_overlap_cutoff', default=0.2, type=float, help="CDS repeats overlap cutoff under which proteins are re-considered as unreliable")
  parser.add_argument('--unreliable_c_score_cutoff', default=0.9, type=float, help="C score cutoff under which proteins deemed unreliable due to repeats overlap are filtered")
  parser.add_argument('--unreliable_prot_coverage_cutoff', default=0.9, type=float, help="Protein coverage cutoff under which proteins deemed unreliable due to repeats overlap are filtered")
  parser.add_argument('--first_command', default=1, type=int, help="Index of first command to be run")
  parser.add_argument('--last_command', default=1, type=int, help="Index of last command to be run")
  args = parser.parse_args()

  main_pipe = Pipe("Annotation_proteins_filtration")
  stages = [blastp,calc_c_score,calc_prot_coverage,calc_repeats_overlap,filter_proteins]

  for i in range(len(stages)):
    if i+1 < args.first_command or i+1 > args.last_command:
      skip = True
    else:
      skip = False
    

