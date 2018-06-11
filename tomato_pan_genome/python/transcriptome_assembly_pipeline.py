"""
TO FIX:
1. Does it still stop after each stage?
3. Why does it proceed to cleanup after BUSCO fails?
5. Integrate STAR and Discasm
"""

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
    r1 = "%s/%s" %(dir_path, r1)
    pos_r2 = r1.replace('_1','_2')
    if os.path.exists(pos_r2):
      PE.append([r1,pos_r2])
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

def transcriptome_assembly_pipeline(data_set_name,sra_accessions,download_target,analysis_target,reference_annotation=None,reference_genome=None, force_overwrite=False, first_command=1, last_command=999):
  """
  
  """
  logging.info("=== TRANSCRIPTOME ASSEMBLY PIPELINE STARTED ===")
  ### 1 - download_data
  logging.info("~~~ STEP 1 - download data ~~~")
  download_command = "%s --gzip --defline-seq '@$sn[_$rn]/$ri' --skip-technical --readids --dumpbase --split-files --clip -O %s %s" % (FASTQ_DUMP_PATH, download_target, sra_accessions)
  if first_command <= 1 and last_command >= 1:
    logging.info("Started downloading data")
    mkdir_overwrite(download_target, overwrite_mode=force_overwrite)
    job_id, exit_status = send_commands_to_queue("fastq-dump_%s" % data_set_name, [download_command],queue_conf)
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
    logging.info("%s PE and %s SE libraries detected." %(len(data_set_PE), len(data_set_SE)))
  
  ### 2 - alignment to reference annotation
  logging.info("~~~ STEP 2 - alignment to reference annotation ~~~")
  if reference_genome:
    topHat_dir = "%s/topHat_alignment" % analysis_target
    genome_index_base = os.path.splitext(reference_genome)[0]
    alignment_commsnds = ["module load python/python-2.7.6",
                "export PATH=/share/apps/bowtie112/bin:$PATH"]
    in_fastq_str = ','.join([l[0] for l in data_set_PE] + data_set_SE) + ' ' + ','.join([l[1] for l in data_set_PE])
    if reference_annotation:
      topHat_command = "%s -p 20 -i 20 -I 5000 -G %s -o %s %s %s" % (TOPHAT_EXEC_PATH, reference_annotation, topHat_dir, genome_index_base, in_fastq_str)
    else:
      topHat_command = "%s -p 20 -i 20 -I 5000 -o %s %s %s" % (TOPHAT_EXEC_PATH, topHat_dir, genome_index_base, in_fastq_str)
    alignment_commsnds.append(topHat_command)
    # sort result
    sort_bam_command = "%s sort %s/accepted_hits.bam %s/accepted_hits.bam.sort" %(SAMTOOLS_EXEC_PATH, topHat_dir, topHat_dir)
    alignment_commsnds.append(sort_bam_command)
    if first_command <= 2 and last_command >= 2:
      logging.info("Started alignment")
      mkdir_overwrite(analysis_target, overwrite_mode=force_overwrite)
      mkdir_overwrite(topHat_dir, overwrite_mode=force_overwrite)
      job_id, exit_status = send_commands_to_queue("topHat_alignment_%s" % data_set_name, alignment_commsnds, queue_conf, n_cpu = 20)
      if exit_status != 0:
        logging.error("Alignment failed. Terminating.")
        sys.exit(1)
    else:
      logging.info("Skipping step...")
  
  ### 3 - transcriptome assembly
  logging.info("~~~ STEP 3 - transcriptome assembly ~~~")
  trinity_assembly_commands = ['bowtie2/bowtie2-2.3.4.1', 'gcc/gcc620', 'htslib/htslib-1.3.2', 'java/java-1.8', 'jellyfish/jellyfish-2.2.7', 'perl/perl-5.20.1-threaded', 'python/anaconda3-5.0.0', 'salomon/salmon-0.9.1', 'samtools/samtools-1.3.1', 'Trinity/Trinity-v2.6.6', 'zlib/zlib129']
  trinity_dir = "%s/trinity_assembly" % analysis_target
  # all paired end
  if data_set_PE and not data_set_SE:
    in_fastq_str = "--left " + ','.join([l[0] for l in data_set_PE]) + " --right " + ','.join([l[1] for l in data_set_PE])  
  # all single end
  elif data_set_SE and not data_set_PE:
    in_fastq_str = "--single " + ','.join(data_set_SE)
  # mixed PE and SE - need to unzip and concat
  else:
    concat_commands = []
    concat_commands.append("zcat -f %s | gzip > %s" %( ' '.join([l[0] for l in data_set_PE] + data_set_SE), download_target + '/concat_R1.fq.gz'))
    concat_commands.append("zcat -f %s | gzip > %s" %( ' '.join([l[1] for l in data_set_PE]), download_target + '/concat_R2.fq.gz'))
    if first_command <= 3 and last_command >= 3:
      logging.info("Handling mixed PE and sE libraries")
      #job_id, exit_status = send_commands_to_queue("%s_concat_reads" % data_set_name,concat_commands,queue_conf)
      exit_status = 0
      if exit_status != 0:
        logging.error("Failed to concatenate reads. Terminating")
        sys.exit(1)
    in_fastq_str = "--left %s --right %s" %(download_target + '/concat_R1.fq.gz', download_target + '/concat_R2.fq.gz')   
  
  if reference_genome:
    sorted_bam = "%s/accepted_hits.bam.sort.bam" % topHat_dir
    assembly_command = "Trinity --seqType fq %s --output %s --CPU 20 --max_memory 50G --genome_guided_max_intron 5000 --genome_guided_bam %s" %(in_fastq_str, trinity_dir, sorted_bam)
    trinity_out = "%s/Trinity-GG.fasta" % trinity_dir
  else:
    assembly_command = "Trinity --seqType fq %s --output %s --CPU 20 --max_memory 50G" %(in_fastq_str, trinity_dir)
    trinity_out = "%s/Trinity.fasta" % trinity_dir
  trinity_assembly_commands.append(assembly_command)
  if first_command <= 3 and last_command >= 3:
    logging.info("Starting assembly")
    os.mkdir(trinity_dir)
    job_id, exit_status = send_commands_to_queue("%s_transcriptome_assembly" % data_set_name, trinity_assembly_commands, queue_conf, n_cpu = 20)
    if exit_status != 0:
      logging.error("Failed in transcriptome assembly. Terminating.")
      sys.exit(1)
  else:
    logging.info("Skipping step...")
    
  ### 4 - BUSCO
  logging.info("~~~ STEP 4 - BUSCO ~~~")
  if first_command <= 4 and last_command >= 4:
    logging.info("Starting BUSCO run")
    busco_commands = ["export PATH=\"/share/apps/augustus/bin:$PATH\"","export PATH=\"/share/apps/augustus/scripts:$PATH\"","export AUGUSTUS_CONFIG_PATH=\"/groups/itay_mayrose/liorglic/software/busco/augustus_config\"", "module load python/python-3.3.0"]
    busco_commands.append("cd %s" % analysis_target)
    busco_commands.append("python %s --in %s --out BUSCO --lineage_path %s --mode transcriptome --cpu 20" %(busco_script_path, trinity_out, busco_lineage_path) )
    job_id, exit_status = send_commands_to_queue("%s_BUSCO" % data_set_name, busco_commands, queue_conf, n_cpu = 20)
    if exit_status != 0:
      logging.error("Failed to run BUSCO. Terminating")
      sys.exit(1)
  else:
    logging.info("Skipping step...")
  
  ### 5 - cleanup
  logging.info("~~~ STEP 5 - Cleanup ~~~")
  if first_command <= 5 and last_command >= 5:
    logging.info("Starting cleanup")
    try:
      # clean reads
      rmtree(download_target)
    except Exception as e:
      logging.warning("Faiiled to clean reads - %s" % str(e))
    try:
      # clean alignment dir
      if reference_genome:
        rmtree(topHat_dir)
    except Exception as e:
      logging.warning("Failed to clean alignment dir - %s" % str(e))
    try:
      # clean assembly dir - delete everything except final output
      for f in os.listdir(trinity_dir):
        if not f.startswith('Trinity') or not f.endswith('.fasta'):
          os.remove("%s/%s" %(trinity_dir,f))
    except Exception as e:
      logging.warning("Failed to clean assembly dir - %s" % str(e))
    try:
      # clean BUSCO dir - remove everything except summary files
      busco_dir = "%s/run_BUSCO/" % analysis_target
      for f in os.listdir(busco_dir):
        if f != 'full_table_BUSCO.tsv' and f != 'short_summary_BUSCO.txt':
          if os.path.isdir("%s/%s" %(busco_dir,f)):
            rmtree("%s/%s" %(busco_dir,f))
          else:
            os.remove("%s/%s" %(busco_dir,f))
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
  TOPHAT_EXEC_PATH = "/share/apps/tophat210/bin/tophat"
  SAMTOOLS_EXEC_PATH = "/share/apps/samtools12/bin/samtools"
  busco_script_path = "/groups/itay_mayrose/liorglic/software/busco/scripts/run_BUSCO.py"
  busco_lineage_path = "/groups/itay_mayrose/liorglic/software/busco/embryophyta_odb9/"

  queue_conf = args.queue_conf
  logging.basicConfig(filename=args.log_file, level=logging.INFO, format='%(asctime)s: %(levelname)s: %(message)s')
  logging.info("Pipeline script command:\npython %s" % ' '.join(sys.argv))

  transcriptome_assembly_pipeline(args.data_set_name, args.sra_accessions, args.download_target, args.analysis_target,args.reference_annotation,args.reference_genome,args.forcie_overwrite,args.first_command,args.last_command)
