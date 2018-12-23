"""
A pipeline for the whole genome annotation procedure, based on MAKER
including the following steps:
1. Annotation lift-over for the official annotation
   This steps results in a gff to be given to next step
2. Full annotation, based on transcript, protein and gene
   model evidence
3. Annotation filtration and refining
4. Compute statistics
5. Cleanup
"""
from __future__ import print_function
import sys
import os
import itertools
import pbs
from shutil import rmtree, copyfile
import argparse
import logging


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

def mkdir_overwrite(path, overwrite_mode):
  try:
    os.makedirs(path)
  except OSError:
    if overwrite_mode:
      rmtree(path)
      os.makedirs(path)
    else:
      raise OSError

def multi_substitute(f_in, f_out, subs_tuple):
  with open(f_in) as f,  open(f_out, 'w') as fo:
    text = f.read()
    subs_text = reduce(lambda a, kv: a.replace(*kv), subs_tuple, text)
    print(subs_text, file=fo, end='')

def genome_annotation_pipeline(genome_name, genome_fasta, out_dir, force_overwrite=False, official_transcripts=None)
  """
  Run genome annotation pipeline
  1. Annotation lift-over for the official annotation
     This steps results in a gff to be given to next step
  2. Full annotation, based on transcript, protein and gene
     model evidence
  3. Annotation filtration and refining
  4. Compute statistics
  5. Cleanup
  """
  logging.info("=== GENOME ANNOTATION PIPELINE STARTED ===")
  mkdir_overwrite(out_dir, force_overwrite)

  ### 1 - Official annotation lift-over
  logging.info("~~~ STEP 1 - Official annotation lift-over ~~~")
  if first_command <= 1 and last_command >= 1:
    if official_transcripts is None:
      logging.error("No official transcripts set provided. You can skip the step using the first_command option")
      sys.exit(1)
    if not os.path.isfile(official_transcripts):
      logging.error("Official transcripts file %s not found." % official_transcripts)
      sys.exit(1)
    # prepare MAKER configurations
    copyfile("%s/maker_bopts.ctl" % TEMPLATES_DIR, out_dir)
    copyfile("%s/maker_exe.ctl" % TEMPLATES_DIR, out_dir)
    liftover_template = "%s/maker_opts_liftover.ctl" % TEMPLATES_DIR
    liftover_maker_conf = out_dir + "/maker_opts.ctl"
    substitutions = ( ('<genome_fasta>', genome_fasta), ('<transcripts_fasta>', official_transcripts) )
    multi_substitute(liftover_template, liftover_maker_conf, substitutions)
    # run MAKER
  else:
    logging.info("Skipping step...")

  ### 2 - Genome annotation
  logging.info("~~~ STEP 2 - Genome annotation ~~~")
  if first_command <= 2 and last_command >= 2:
    pass
  else:
    logging.info("Skipping step...")

  ### 3 - Annotation filtration and refining
  logging.info("~~~ STEP 3 - Annotation filtration and refining ~~~")
  if first_command <= 3 and last_command >= 4:
    pass
  else:
    logging.info("Skipping step...")

  ### 4 - Compute statistics
  logging.info("~~~ STEP 1 - Compute statistics ~~~")
  if first_command <= 4 and last_command >= 4:
    pass
  else:
    logging.info("Skipping step...")

  ### 5 - Cleanup
  logging.info("~~~ STEP 5 - Cleanup ~~~")
  if first_command <= 5 and last_command >= 5:
    pass
  else:
    logging.info("Skipping step...")

  logging.info("=== GENOME ANNOTATION PIPELINE COMPLETED SUCCESSFULLY ===")


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('genome_name', help="Identifier for the genome to be annotated")
  parser.add_argument('genome_fasta', help="Path to genome fasta to be annotated", action=FullPaths)
  parser.add_argument('out_dir', help="Path to output directory", action=FullPaths)
  parser.add_argument('log_file', action=FullPaths, help="Path to run log file")
  parser.add_argument('--official_transcripts_set', default=None, action=FullPaths, help="Path to fasta file with transcripts derived from official annotation")
  parser.add_argument('--full_transcripts_set', default=None, action=FullPaths, help="Path to fasta file with full transcripts set")
  parser.add_argument('--full_proteins_set', default=None, action=FullPaths, help="Path to fasta file with full proteins set")
  parser.add_argument('-f', '--force_overwrite', action='store_true', default=False)
  parser.add_argument('--first_command', default=1, type=int)
  parser.add_argument('--last_command', default=999, type=int)
  args = parser.parse_args()

  BUSCO_SCRIPT_PATH = "/groups/itay_mayrose/liorglic/software/busco/scripts/run_BUSCO.py"
  BUSCO_LINEAGE_PATH = "/groups/itay_mayrose/liorglic/software/busco/embryophyta_odb9/"
  TEMPLATES_DIR = ""	# for MAKER configuration templates

  # set logger, including redirect STDOUT and STDERR to log
  logging.basicConfig(filename=args.log_file, level=logging.INFO,
                      format='%(asctime)s: %(name)s:%(levelname)s: %(message)s')
  stdout_logger = logging.getLogger('STDOUT')
  sl = StreamToLogger(stdout_logger, logging.INFO)
  sys.stdout = sl
  stderr_logger = logging.getLogger('STDERR')
  sl = StreamToLogger(stderr_logger, logging.ERROR)
  sys.stderr = sl

  logging.info("Pipeline script command:\npython %s" % ' '.join(sys.argv))

  # run pipeline
  genome_annotation_pipeline()
