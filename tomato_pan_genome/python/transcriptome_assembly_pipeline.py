import sys
import os

sys.path.append("../../queue_utilities/")
from queueUtils import send_commands_to_queue
from shutil import rmtree
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


def prep_libs_lists(dir_path):
  """
  Parses file names in dir and returns a tuple with two elements:
  0. a list of lists of PE files
  1. a list of SE files
  """
  PE = []
  SE = []
  dir_files = os.listdir(dir_path)
  r1_files = [f for f in dir_files if
              (f.endswith('fastq.gz') or f.endswith('fq.gz') or f.endswith('fastq') or f.endswith('fq'))
              and '_1' in f]
  for r1 in r1_files:
    r1 = "%s/%s" % (dir_path, r1)
    pos_r2 = r1.replace('_1', '_2')
    if os.path.exists(pos_r2):
      PE.append([r1, pos_r2])
    else:
      SE.append(r1)
  return PE, SE


def mkdir_overwrite(path, overwrite_mode):
  try:
    os.makedirs(path)
  except OSError:
    if overwrite_mode:
      rmtree(path)
      os.makedirs(path)
    else:
      raise OSError

def transcriptome_assembly_pipeline(data_set_name, sra_accessions, download_target, analysis_target,
                                    reference_annotation=None, reference_genome=None, force_overwrite=False,
                                    first_command=1, last_command=999):
  """

  """
  logging.info("=== TRANSCRIPTOME ASSEMBLY PIPELINE STARTED ===")
  ### 1 - download_data
  logging.info("~~~ STEP 1 - download data ~~~")
  download_command = "%s --gzip --defline-seq '@$sn[_$rn]/$ri' --skip-technical --readids --dumpbase --split-files --clip -O %s %s" % (
    FASTQ_DUMP_PATH, download_target, sra_accessions)
  if first_command <= 1 and last_command >= 1:
    logging.info("Started downloading data")
    mkdir_overwrite(download_target, overwrite_mode=force_overwrite)
    job_id, exit_status = send_commands_to_queue("fastq-dump_%s" % data_set_name, [download_command], queue_conf)
    if exit_status != 0:
      logging.error("Download failed. Terminating")
      sys.exit(1)
  else:
    logging.info("Skipping step...")

  # parse libs
  data_set_PE, data_set_SE = prep_libs_lists(download_target)
  if not data_set_PE and not data_set_SE:
    logging.error("No files found for data set %s. Terminating." % data_set_name)
    sys.exit(1)
  else:
    logging.info("%s PE and %s SE libraries detected." % (len(data_set_PE), len(data_set_SE)))

  ### 2 - alignment to reference annotation
  logging.info("~~~ STEP 2 - alignment to reference genome ~~~")
  if reference_genome:
    STAR_dir = "%s/STAR_alignment" % analysis_target
    alignment_commands = ["module load gcc/gcc620 perl/perl518 STAR/STAR-2.6.0a"]
    # SE alignment
    if data_set_SE:
      in_fastq_str = ','.join(data_set_SE)
      alignment_commands.append(
        "STAR --genomeDir %s --readFilesIn %s --twopassMode Basic --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --runThreadN 5 --outSAMstrandField intronMotif --outFileNamePrefix %s/SE_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --readFilesCommand zcat" % (
          reference_genome, in_fastq_str, STAR_dir))
    se_mapped_bam = "%s/SE_Aligned.sortedByCoord.out.bam" % STAR_dir
    se_unmapped_out = "%s/SE_Unmapped.out.mate1" % STAR_dir
    se_unmapped_fq = "%s/SE_Unmapped_R1.fq" % STAR_dir
    alignment_commands.append(r"sed 's/\(@SOLEXA.*\)/\1\/1/' %s > %s" % (se_unmapped_out, se_unmapped_fq))
    # PE alignment
    if data_set_PE:
      in_fastq_str = ','.join([l[0] for l in data_set_PE]) + ' ' + ','.join([l[1] for l in data_set_PE])
      alignment_commands.append(
        "STAR --genomeDir %s --readFilesIn %s --twopassMode Basic --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --runThreadN 5 --outSAMstrandField intronMotif --outFileNamePrefix %s/PE_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --readFilesCommand zcat" % (
          reference_genome, in_fastq_str, STAR_dir))
      pe_mapped_bam = "%s/PE_Aligned.sortedByCoord.out.bam" % STAR_dir
      pe_unmapped_out1 = "%s/PE_Unmapped.out.mate1" % STAR_dir
      pe_unmapped_out2 = "%s/PE_Unmapped.out.mate2" % STAR_dir
      pe_unmapped_fq1 = "%s/PE_Unmapped_R1.fq" % STAR_dir
      pe_unmapped_fq2 = "%s/PE_Unmapped_R2.fq" % STAR_dir
      alignment_commands.append(r"sed 's/\t.*/\/1/' %s > %s" % (pe_unmapped_out1, pe_unmapped_fq1))
      alignment_commands.append(r"sed 's/\t.*/\/2/' %s > %s" % (pe_unmapped_out2, pe_unmapped_fq2))
    # define final outputs and combine SE and PE results, if needed
    if data_set_SE and not data_set_PE:
      final_aligned_bam = se_mapped_bam
      final_unmapped_fq1 = se_unmapped_fq
      final_unmapped_fq2 = None
    elif data_set_PE and not data_set_SE:
      final_aligned_bam = pe_mapped_bam
      final_unmapped_fq1 = pe_unmapped_fq1
      final_unmapped_fq2 = pe_unmapped_fq2
    elif data_set_PE and data_set_SE:
      final_aligned_bam = "%s/combined_Aligned.sortedByCoord.out.bam" % STAR_dir
      final_unmapped_fq1 = "%s/combined_Unmapped_R1.fq" % STAR_dir
      final_unmapped_fq2 = pe_unmapped_fq2
      alignment_commands.append(
        "%s merge %s %s %s" % (SAMTOOLS_EXEC_PATH, final_aligned_bam, pe_mapped_bam, se_mapped_bam))
      alignment_commands.append("cat %s %s > %s" % (pe_unmapped_fq1, se_unmapped_fq, final_unmapped_fq1))

    if first_command <= 2 and last_command >= 2:
      logging.info("Started alignment")
      mkdir_overwrite(analysis_target, overwrite_mode=force_overwrite)
      mkdir_overwrite(STAR_dir, overwrite_mode=force_overwrite)
      job_id, exit_status = send_commands_to_queue("STAR_alignment_%s" % data_set_name, alignment_commands,
                                                   queue_conf, n_cpu=5)
      if exit_status != 0:
        logging.error("Alignment failed. Terminating.")
        sys.exit(1)
    else:
      logging.info("Skipping step...")
    # ensure final results exist
    if last_command >= 2:
      missing_outputs = False
      final_out_files = (final_aligned_bam, final_unmapped_fq1, final_unmapped_fq2)
      for f in final_out_files:
        if f and not os.path.exists(f):
          logging.error("Expected output file %s not found. Terminating." % f)
          missing_outputs = True
      if missing_outputs:
        sys.exit(1)

  ### 3 - transcriptome assembly
  logging.info("~~~ STEP 3 - transcriptome assembly ~~~")
  trinity_assembly_commands = ['module unload R/R301 gcc/gcc480 python/python-3.3.0',
                               'module load salomon/salmon-0.9.1 perl/perl-5.20.1-threaded Trinity/Trinity-v2.6.6 zlib/zlib129 java/java-1.8 gcc/gcc620 samtools/samtools-1.3.1 htslib/htslib-1.3.2 jellyfish/jellyfish-2.2.7 python/anaconda3-5.0.0',
                               'module load perl/perl-5.20.1-threaded zlib/zlib129 jellyfish/jellyfish-2.2.7 gcc/gcc620 samtools/samtools-1.3.1 htslib/htslib-1.3.2 salomon/salmon-0.9.1',
                               'module load Trinity/Trinity-v2.6.6', 'module load java/java-1.8',
                               'module load jellyfish', 'module load bowtie2/bowtie2-2.3.4.1',
                               'export PATH=/groups/itay_mayrose/liorglic/Salmon-latest_linux_x86_64/bin:/groups/itay_mayrose/liorglic/samtools-1.7/bin:/share/apps/Jellyfish-2.2.7/bin:$PATH']
  trinity_dir_mapped = "%s/trinity_assembly_mapped" % analysis_target
  trinity_dir_unmapped = "%s/trinity_assembly_unmapped" % analysis_target
  # all paired end
  if data_set_PE and not data_set_SE:
    in_fastq_str = "--left " + ','.join([l[0] for l in data_set_PE]) + " --right " + ','.join(
      [l[1] for l in data_set_PE])
  # all single end
  elif data_set_SE and not data_set_PE:
    in_fastq_str = "--single " + ','.join(data_set_SE)
  # mixed PE and SE - need to unzip and concat
  else:
    concat_commands = []
    concat_commands.append("zcat -f %s | gzip > %s" % (
      ' '.join([l[0] for l in data_set_PE] + data_set_SE), download_target + '/concat_R1.fq.gz'))
    concat_commands.append(
      "zcat -f %s | gzip > %s" % (' '.join([l[1] for l in data_set_PE]), download_target + '/concat_R2.fq.gz'))
    if first_command <= 3 and last_command >= 3:
      logging.info("Handling mixed PE and sE libraries")
      job_id, exit_status = send_commands_to_queue("%s_concat_reads" % data_set_name, concat_commands, queue_conf)
      exit_status = 0
      if exit_status != 0:
        logging.error("Failed to concatenate reads. Terminating")
        sys.exit(1)
    in_fastq_str = "--left %s --right %s" % (
      download_target + '/concat_R1.fq.gz', download_target + '/concat_R2.fq.gz')

  if reference_genome:
    # assembly of mapped reads, using reference
    trinity_assembly_commands.append(
      "Trinity --seqType fq %s --output %s --CPU 20 --max_memory 50G --genome_guided_max_intron 5000 --genome_guided_bam %s" % (
        in_fastq_str, trinity_dir_mapped, final_aligned_bam))
    # assembly of unmapped reads
    if data_set_PE:
      in_fastq_str = "--left %s --right %s" % (final_unmapped_fq1, final_unmapped_fq2)
    else:
      in_fastq_str = "--left %s" % (final_unmapped_fq1)
    trinity_assembly_commands.append(
      "Trinity --seqType fq %s --output %s --CPU 20 --max_memory 50G" % (in_fastq_str, trinity_dir_unmapped))
    # concat transcripts from mapped and unmapped runs
    trinity_out = "%s/Trinity-combined.fasta" % analysis_target
    trinity_assembly_commands.append(
      "cat %s/Trinity-GG.fasta %s/Trinity.fasta > %s" % (trinity_dir_mapped, trinity_dir_unmapped, trinity_out))
  else:
    trinity_assembly_commands.append(
      "Trinity --seqType fq %s --output %s --CPU 20 --max_memory 50G" % (in_fastq_str, trinity_dir_unmapped))
    trinity_out = "%s/Trinity.fasta" % trinity_dir_unmapped
  if first_command <= 3 and last_command >= 3:
    logging.info("Starting assembly")
    if data_set_PE:
      os.mkdir(trinity_dir_mapped)
    os.mkdir(trinity_dir_unmapped)
    job_id, exit_status = send_commands_to_queue("%s_transcriptome_assembly" % data_set_name,
                                                 trinity_assembly_commands, queue_conf, n_cpu=20)
    if exit_status != 0:
      logging.error("Failed in transcriptome assembly. Terminating.")
      sys.exit(1)
  else:
    logging.info("Skipping step...")
  # ensure output was created
  if last_command >= 3 and not os.path.exists(trinity_out):
    logging.error("Expected output file %s not found. Terminating." % trinity_out)
    sys.exit(1)

  ### 4 - BUSCO
  logging.info("~~~ STEP 4 - BUSCO ~~~")
  if first_command <= 4 and last_command >= 4:
    logging.info("Starting BUSCO run")
    busco_commands = ["export PATH=\"/share/apps/augustus/bin:$PATH\"",
                      "export PATH=\"/share/apps/augustus/scripts:$PATH\"",
                      "export AUGUSTUS_CONFIG_PATH=\"/groups/itay_mayrose/liorglic/software/busco/augustus_config\"",
                      "module load python/python-3.3.0"]
    busco_commands.append("cd %s" % analysis_target)
    busco_commands.append("python %s --in %s --out BUSCO --lineage_path %s --mode transcriptome --cpu 20" % (
      busco_script_path, trinity_out, busco_lineage_path))
    job_id, exit_status = send_commands_to_queue("%s_BUSCO" % data_set_name, busco_commands, queue_conf, n_cpu=20)
    if exit_status != 0:
      logging.error("Failed to run BUSCO. Terminating")
      sys.exit(1)
  else:
    logging.info("Skipping step...")
  # ensure final output exists
  busco_out = "%s/run_BUSCO/short_summary_BUSCO.txt" % analysis_target
  if last_command >= 4 and not os.path.exists(busco_out):
    logging.error("Expected output %s not found. Terminating." % busco_out)
    sys.exit(1)

  ### 5 - Collect stats
  logging.info("~~~ STEP 5 - Collect stats ~~~")
  if first_command <= 5 and last_command >= 5:
    logging.info("Starting to collect stats")
    stats_out = "%s/transcriptome_stats.tsv" % analysis_target
    collect_stats_commands = ["module load python/python-2.7.6", "python %s/get_genome_stats.py %s > %s" %(os.path.dirname(os.path.realpath(__file__)), trinity_out, stats_out)]
    job_id, exit_status = send_commands_to_queue("%s_collect_stats" % data_set_name, collect_stats_commands, queue_conf)
    if exit_status != 0:
      logging.error("Failed to collect stats. Terminating")
      sys.exit(1)
  else:
    logging.info("Skipping step...")
  # ensure final output exists
  if last_command >= 5 and not os.path.exists(stats_out):
    logging.error("Expected output %s not found. Terminating." % stats_out)
    sys.exit(1)

  ### 6 - cleanup
  logging.info("~~~ STEP 6 - Cleanup ~~~")
  if first_command <= 6 and last_command >= 6:
    logging.info("Starting cleanup")
    try:
      # clean reads
      rmtree(download_target)
    except Exception as e:
      logging.warning("Faiiled to clean reads - %s" % str(e))
    try:
      # clean alignment dir
      if reference_genome:
        rmtree(STAR_dir)
    except Exception as e:
      logging.warning("Failed to clean alignment dir - %s" % str(e))
    try:
      # clear tmp dir
      rmtree("%s/tmp" % analysis_target)
    except Exception as e:
      logging.warning("Failed to clean tmp dir - %s" % str(e))
    try:
      # clean assembly dir - delete everything except final output
      for f in os.listdir(trinity_dir_unmapped):
        if not f.startswith('Trinity') or not f.endswith('.fasta'):
          os.remove("%s/%s" % (trinity_dir_unmapped, f))
      if reference_genome:
        for f in os.listdir(trinity_dir_mapped):
          if not f.startswith('Trinity') or not f.endswith('.fasta'):
            os.remove("%s/%s" % (trinity_dir_mapped, f))
    except Exception as e:
      logging.warning("Failed to clean assembly dir - %s" % str(e))
    try:
      # clean BUSCO dir - remove everything except summary files
      busco_dir = "%s/run_BUSCO/" % analysis_target
      for f in os.listdir(busco_dir):
        if f != 'full_table_BUSCO.tsv' and f != 'short_summary_BUSCO.txt':
          if os.path.isdir("%s/%s" % (busco_dir, f)):
            rmtree("%s/%s" % (busco_dir, f))
          else:
            os.remove("%s/%s" % (busco_dir, f))
    except Exception as e:
      logging.warning("Failed to clean BUSCO dir - %s" % str(e))
  else:
    logging.info("Skipping step...")

  logging.info("i=== TRANSCRIPTOME ASSEMBLY PIPELINE COMPLETED SUCCESSFULLY ===")


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('queue_conf', action=FullPaths)
  parser.add_argument('data_set_name')
  parser.add_argument('sra_accessions')
  parser.add_argument('download_target', action=FullPaths)
  parser.add_argument('analysis_target', action=FullPaths)
  parser.add_argument('log_file', action=FullPaths)
  parser.add_argument('-a', '--reference_annotation', default=None, action=FullPaths)
  parser.add_argument('-g', '--reference_genome', default=None, action=FullPaths)
  parser.add_argument('-f', '--forcie_overwrite', action='store_true', default=False)
  parser.add_argument('--first_command', default=1, type=int)
  parser.add_argument('--last_command', default=999, type=int)
  args = parser.parse_args()

  FASTQ_DUMP_PATH = "/groups/itay_mayrose/liorglic/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump"
  SAMTOOLS_EXEC_PATH = "/share/apps/samtools12/bin/samtools"
  busco_script_path = "/groups/itay_mayrose/liorglic/software/busco/scripts/run_BUSCO.py"
  busco_lineage_path = "/groups/itay_mayrose/liorglic/software/busco/embryophyta_odb9/"

  queue_conf = args.queue_conf
  # set logger, including recirect STDOUT and STDERR to log
  logging.basicConfig(filename=args.log_file, level=logging.INFO,
                      format='%(asctime)s: %(name)s:%(levelname)s: %(message)s')
  stdout_logger = logging.getLogger('STDOUT')
  sl = StreamToLogger(stdout_logger, logging.INFO)
  sys.stdout = sl
  stderr_logger = logging.getLogger('STDERR')
  sl = StreamToLogger(stderr_logger, logging.ERROR)
  sys.stderr = sl

  logging.info("Pipeline script command:\npython %s" % ' '.join(sys.argv))
  transcriptome_assembly_pipeline(args.data_set_name, args.sra_accessions, args.download_target, args.analysis_target,
                                  args.reference_annotation, args.reference_genome, args.forcie_overwrite,
                                  args.first_command, args.last_command)

