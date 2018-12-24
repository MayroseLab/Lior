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
from __future__ import print_function, division
import sys
import os
import itertools
from pbs.job import Job
from shutil import rmtree, copyfile, move
import argparse
import logging
import subprocess

### GLOBALS
MAKER_CPUS = 160
MAKER_NODES = 8
MAKER_MEM = "20g"

### UTILITY CLASSES
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

### FUNCTIONS
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
    with open(f_in) as f, open(f_out, 'w') as fo:
        text = f.read()
        subs_text = reduce(lambda a, kv: a.replace(*kv), subs_tuple, text)
        print(subs_text, file=fo, end='')

def check_maker_run_complete(run_log):
    """
    Takes a maker run log (stderr) and checks if it ends with
    the line "Maker is now finished!!!"
    Return True/False
    """
    p = subprocess.Popen(['tail', run_log], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    line = p.stdout.readline()
    while line:
        if line.startswith("Maker is now finished!!!"):
            return True
    return False

def genome_annotation_pipeline(genome_name, genome_fasta, out_dir, config_templates, force_overwrite=False,
                               official_transcripts=None, annotation_transcripts='',
                               annotation_proteins='', external_gff='', dryrun=False, first_command=1, last_command=999):
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
    genome_fasta_base = os.path.splitext(os.path.basename(genome_fasta))[0]

    ### 1 - Official annotation lift-over
    logging.info("~~~ STEP 1 - Official annotation lift-over ~~~")
    if first_command <= 1 and last_command >= 1:
        if official_transcripts is None:
            logging.error(
                "No official transcripts set provided. You can skip the "
                "step using the first_command option")
            sys.exit(1)
        if not os.path.isfile(official_transcripts):
            logging.error("Official transcripts file %s not found." % official_transcripts)
            sys.exit(1)
        # prepare MAKER configurations
        copyfile("%s/maker_bopts.ctl" % config_templates, out_dir)
        copyfile("%s/maker_exe.ctl" % config_templates, out_dir)
        liftover_template = "%s/maker_opts_liftover.ctl" % config_templates
        liftover_maker_conf = out_dir + "/maker_opts.ctl"
        substitutions = (
        ('<genome_fasta>', genome_fasta), ('<transcripts_fasta>', official_transcripts))
        multi_substitute(liftover_template, liftover_maker_conf, substitutions)
        # run MAKER
        liftover_base = genome_fasta_base + "_liftover"
        commands = ["module load miniconda/miniconda2-4.5.4-MakerMPI",
                    "cd %s" % out_dir,
                    "mpirun -n %s -env I_MPI_FABRICS tcp maker -base %s > maker_liftover.out 2> "
                    "maker_liftover.err" % (MAKER_CPUS, liftover_base)]
        job_name = "%s_maker_liftover" % genome_name
        maker_liftover_job = Job(job_name, commands=commands,
                                 nodes=8, ppn=int(MAKER_CPUS/MAKER_NODES), pmem=MAKER_MEM)
        # write script to file
        maker_liftover_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue (block)
            logging.info("Running MAKER (liftover)...")
            maker_liftover_job.submit_block()
            # check that run completed successfully
            if not check_maker_run_complete("%s/maker_liftover.err" % out_dir):
                logging.error("MAKER run failed (see error in %s/maker_liftover.err). "
                              "Terminating pipeline" % out_dir)
                sys.exit(1)
        # make GFF
        liftover_out_dir = "%s/%s.maker.output/" %(out_dir, liftover_base)
        liftover_master_index = "%s/%s_master_datastore_index.log" %(liftover_out_dir,
                                                                     liftover_base)
        commands = ["module load miniconda/miniconda2-4.5.4-MakerMPI",
                    "cd %s" % liftover_out_dir,
                    "gff3_merge -d %s -g -n "
                    "> %s/liftover_merge_gff.out 2> %s/liftover_merge_gff.err"
                    %(liftover_master_index, out_dir, out_dir)]
        job_name = "%s_liftover_merge_gff" % genome_name
        liftover_make_gff_job = Job(job_name, commands=commands)
        # write script to file
        liftover_make_gff_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue (block)
            logging.info("Creating merged GFF...")
            liftover_make_gff_job.submit_block()
    else:
        logging.info("Skipping step...")
    liftover_gff = "%s/%s.all.gff" %(liftover_out_dir, liftover_base)

    ### 2 - Genome annotation
    logging.info("~~~ STEP 2 - Genome annotation ~~~")
    if first_command <= 2 and last_command >= 2:
        # prepare MAKER configurations
        copyfile("%s/maker_bopts.ctl" % config_templates, out_dir)
        copyfile("%s/maker_exe.ctl" % config_templates, out_dir)
        annotation_template = "%s/maker_opts.ctl" % config_templates
        annotation_maker_conf = out_dir + "/maker_opts.ctl"
        substitutions = [
            ('<genome_fasta>', genome_fasta),
            ('<transcripts_fasta>', annotation_transcripts),
            ('<proteins_fasta>', annotation_proteins),
        ]
        gff_input = ''
        if os.path.exists(liftover_gff):
            gff_input += liftover_gff
        if external_gff:
            gff_input += "," + external_gff
        substitutions.append(('<gene_models_gff>', gff_input))
        multi_substitute(annotation_template, annotation_maker_conf, substitutions)
        # run MAKER
        commands = ["module load miniconda/miniconda2-4.5.4-MakerMPI",
                    "cd %s" % out_dir,
                    "mpirun -n %s -env I_MPI_FABRICS tcp maker > maker_annotation.out 2> "
                    "maker_annotation.err" % MAKER_CPUS]
        job_name = "%s_maker_annotation" % genome_name
        maker_annotation_job = Job(job_name, commands=commands,
                                 nodes=8, ppn=int(MAKER_CPUS/MAKER_NODES), pmem=MAKER_MEM)
        # write script to file
        maker_annotation_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue (block)
            logging.info("Running MAKER (annotation)...")
            maker_annotation_job.submit_block()
            # check that run completed successfully
            if not check_maker_run_complete("%s/maker_annotation.err" % out_dir):
                logging.error("MAKER run failed (see error in %s/maker_annotation.err). "
                              "Terminating pipeline" % out_dir)
                sys.exit(1)
        # make merged gff and fasta
        annotation_out_dir = "%s/%s.maker.output/" %(out_dir, genome_fasta_base)
        annotation_master_index = "%s/%s_master_datastore_index.log" %(annotation_out_dir,
                                                                     genome_fasta_base)
        commands = ["module load miniconda/miniconda2-4.5.4-MakerMPI",
                    "cd %s" % annotation_out_dir,
                    "gff3_merge -d %s -g -n "
                    "> %s/annotation_merge_gff.out 2> %s/annotation_merge_gff.err"
                    %(annotation_master_index, out_dir, out_dir),
                    "fasta_merge -d %s" % annotation_master_index]
        job_name = "%s_annotation_merge" % genome_name
        annotation_merge_job = Job(job_name, commands=commands)
        # write script to file
        annotation_merge_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue (block)
            logging.info("Running gff and fasta merge...")
            logging.info("Creating merged GFF and FASTA...")
            annotation_merge_job.submit_block()
    else:
        logging.info("Skipping step...")
    annotation_raw_gff = "%s/%s.all.gff" %(annotation_out_dir, genome_fasta_base)
    annotation_raw_fasta = "%s/%s.all.maker.proteins.fasta" %(annotation_out_dir, genome_fasta_base)

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
    parser.add_argument('genome_fasta', help="Path to genome fasta to be annotated",
                        action=FullPaths)
    parser.add_argument('out_dir', help="Path to output directory", action=FullPaths)
    parser.add_argument('log_file', action=FullPaths, help="Path to run log file")
    parser.add_argument('config_templates', action=FullPaths, help="Path to configurations templates dir")
    parser.add_argument('--official_transcripts_set', default=None, action=FullPaths,
                        help="Path to fasta file with transcripts derived from official annotation")
    parser.add_argument('--full_transcripts_set', default=None, action=FullPaths,
                        help="Path to fasta file with full transcripts set")
    parser.add_argument('--full_proteins_set', default=None, action=FullPaths,
                        help="Path to fasta file with full proteins set")
    parser.add_argument('--external_gff', default=None, action=FullPaths,
                        help="Path to an external gff with gene models to be used")
    parser.add_argument('-f', '--force_overwrite', action='store_true', default=False)
    parser.add_argument('--first_command', default=1, type=int, help="First command index ("
                                                                     "1-based)")
    parser.add_argument('--last_command', default=999, type=int, help="Last command index (1-based")
    parser.add_argument('--dryrun', action="store_true", default=False, help="Do not send any "
                                                                             "commands to queue")
    args = parser.parse_args()

    BUSCO_SCRIPT_PATH = "/groups/itay_mayrose/liorglic/software/busco/scripts/run_BUSCO.py"
    BUSCO_LINEAGE_PATH = "/groups/itay_mayrose/liorglic/software/busco/embryophyta_odb9/"

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
    genome_annotation_pipeline(args.genome_name, args.genome_fasta, args.out_dir, args.config_templates, force_overwrite=args.force_overwrite,
                               official_transcripts=args.official_transcripts_set, annotation_transcripts=args.full_transcripts_set,
                               annotation_proteins=args.full_proteins_set, external_gff=args.external_gff, dryrun=args.dryrun,
                               first_command=args.first_command, last_command=args.last_command)
