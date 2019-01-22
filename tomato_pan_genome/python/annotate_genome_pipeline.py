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
from shutil import rmtree, copy, move
import argparse
import logging
import subprocess

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
        line = p.stdout.readline()
    return False

def genome_annotation_pipeline(**kwargs):
    """
    Run genome annotation pipeline
    1. Annotation lift-over for the official annotation
       This steps results in a gff to be given to next step
    2. Full annotation, based on transcript, protein and gene
       model evidence
    3. Annotation QA, filtration and refining
    4. Compute statistics
    5. Cleanup
    """
    # fetch args
    for name, value in kwargs.items():
        if type(value) == str:
            exec("%s = \"%s\"" % (name, value))
        else:
            exec("%s = %s" % (name, value))

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
        copy("%s/maker_bopts.ctl" % config_templates, out_dir)
        copy("%s/maker_exe.ctl" % config_templates, out_dir)
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
                    "maker_liftover.err" % (cpus, liftover_base)]
        job_name = "%s_maker_liftover" % genome_name
        maker_liftover_job = Job(job_name, command=commands,
                                 nodes=nodes, ppn=int(cpus/nodes), pmem=maker_ram)
        # write script to file
        maker_liftover_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue (block)
            logging.info("Running MAKER (liftover)...")
            maker_liftover_job.submit_block()
            # check that run completed successfully
            logging.info("Verifying successfull MAKER run...")
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
        liftover_make_gff_job = Job(job_name, command=commands)
        # write script to file
        liftover_make_gff_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue (block)
            logging.info("Creating merged GFF...")
            liftover_make_gff_job.submit_block()
    else:
        logging.info("Skipping step...")
    liftover_gff = "%s/%s.all.gff" %(liftover_out_dir, liftover_base)
    move(liftover_maker_conf, "%s//maker_opts_liftover.ctl" % out_dir)

    ### 2 - Genome annotation
    logging.info("~~~ STEP 2 - Genome annotation ~~~")
    if first_command <= 2 and last_command >= 2:
        # prepare MAKER configurations
        copy("%s/maker_bopts.ctl" % config_templates, out_dir)
        copy("%s/maker_exe.ctl" % config_templates, out_dir)
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
                    "maker_annotation.err" % cpus]
        job_name = "%s_maker_annotation" % genome_name
        maker_annotation_job = Job(job_name, command=commands,
                                 nodes=nodes, ppn=int(cpus/nodes), pmem=maker_ram)
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
                    "gff3_merge -d %s -n "
                    "> %s/annotation_merge_gff.out 2> %s/annotation_merge_gff.err"
                    %(annotation_master_index, out_dir, out_dir),
                    "fasta_merge -d %s" % annotation_master_index]
        job_name = "%s_annotation_merge" % genome_name
        annotation_merge_job = Job(job_name, command=commands)
        # write script to file
        annotation_merge_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue (block)
            logging.info("Creating merged GFF and FASTA...")
            annotation_merge_job.submit_block()
    else:
        logging.info("Skipping step...")
    annotation_raw_gff = "%s/%s.all.gff" %(annotation_out_dir, genome_fasta_base)
    annotation_raw_fasta = "%s/%s.all.maker.proteins.fasta" %(annotation_out_dir, genome_fasta_base)

    ### 3 - Annotation QA, filtration and refining
    logging.info("~~~ STEP 3 - Annotation QA, filtration and refining ~~~")
    if first_command <= 3 and last_command >= 4:
        pass
        ## BUSCO
        busco_cpus = min(10,max(int(cpus/10),1))
        job_name = "%s_BUSCO" % genome_name
        busco_commands = ["export AUGUSTUS_CONFIG_PATH=\"%s\"" % augustus_conf,
                          "export PATH=\"%s:$PATH\"" % augustus_bin,
                          "cd %s" % out_dir,
                          "python %s --in %s --out BUSCO --lineage_path %s --mode proteins --cpu %s >%s/%s.out 2>%s/%s.err"
                          %(busco_exe, annotation_raw_fasta, busco_set, busco_cpus, out_dir, job_name, out_dir, job_name)]
        busco_job = Job(job_name, command=busco_commands, ppn=busco_cpus, nodes=1)
        # write script to file
        busco_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue
            logging.info("Running BUSCO...")
            busco_job.submit_block()
        busco_full_report = "%s/run_BUSCO/full_table_BUSCO.tsv" % out_dir
        ## BLAST - find sequence similarities
        blast_cpus = min(10,max(int(cpus/10),1)) 
        blast_out = "%s/%s_vs_%s.blast.tsv" %(out_dir, genome_name, os.path.basename(blast_qa_db))
        job_name = "%s_blast_vs_proteins_db" % genome_name
        blast_commands = ["module load blast/blast-2.7.1",
                          "blastp -query %s -db %s -evalue 1e-5 -outfmt \"6 qseqid sseqid pident qlen length mismatch gapopen evalue bitscore qcovs\" -num_threads %s -out %s -max_target_seqs 1 >%s/%s.out 2>%s/%s.err" %(annotation_raw_fasta, blast_qa_db, blast_cpus, blast_out, out_dir, job_name, out_dir, job_name)]
        blast_job = Job(job_name, command=blast_commands, ppn=blast_cpus, nodes=1)
        # write script to file
        blast_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue
            logging.info("Running BLAST vs. QA DB...")
            blast_job.submit_block()
        ## InterProScan (ips)
        ips_cpus = min(20,max(int(cpus*0.8),1))
        ips_out_pref = "%s/%s.interProScan" %(out_dir, genome_name)
        job_name = "%s_interProScan" % genome_name
        ips_commands = ["module load interproscan/5.32-71",
                        "module load java/jdk1.8.25",
                        "interproscan.sh -i %s -b %s -t p -dp -pa -appl Pfam,ProDom,SuperFamily,PIRSF --goterms --iprlookup -f tsv,gff3 -T %s >%s/%s.out 2>%s/%s.err" %(annotation_raw_fasta, ips_out_pref, out_dir, out_dir, job_name, out_dir, job_name)]
        ips_job = Job(job_name, command=ips_commands, ppn=ips_cpus, nodes=1)
        # write script to file
        ips_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue
            logging.info("Running InterProScan...")
            ips_job.submit_block()
        ips_out_tsv = ips_out_pref + ".tsv"
        ## Ensure BUSCO, Blast and InterproScan completed successfully
        if not dryrun and not (os.path.isfile(busco_full_report) and os.path.isfile(blast_out) and os.path.isfile(ips_out_tsv)):
            logging.error("One of the QA results is missing. Terminating.")
            sys.exit(1)
        ## create QA report
        job_name = "%s_QA_report" % genome_name
        out_qa_report = "%s/%s.qa_report.tsv" %(out_dir, genome_name)
        qa_report_commands = ["module load python/python-anaconda2.7",
                        "python %s %s %s --busco_result %s --blast_result %s --ips_result %s" %(qa_report_script, annotation_raw_gff, out_qa_report, busco_full_report, blast_out, ips_out_tsv)]
        qa_report_job = Job(job_name, command=qa_report_commands, ppn=1, nodes=1)
        # write script to file
        qa_report_job.script("%s/%s.q" %(out_dir, job_name))
        if not dryrun:
            # send commands to queue
            logging.info("Creating QA report...")
            qa_report_job.submit_block()
        ## Filter annotation
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
    parser.add_argument('--genome_name', required=True, help="Identifier for the genome to be annotated")
    parser.add_argument('--genome_fasta', required=True, help="Path to genome fasta to be annotated",
                        action=FullPaths)
    parser.add_argument('--out_dir', required=True, help="Path to output directory", action=FullPaths)
    parser.add_argument('--log_file', required=True, action=FullPaths, help="Path to run log file")
    parser.add_argument('--config_templates', required=True, action=FullPaths, help="Path to configurations templates dir")
    parser.add_argument('--cpus', type=int, default=1, help="Number of CPUs to use")
    parser.add_argument('--nodes', type=int, default=1, help="Number compute nodes to spread CPUs on")
    parser.add_argument('--maker_ram', default='20g', help="Max RAM to be used by MAKER")
    parser.add_argument('--official_transcripts', default=None, action=FullPaths,
                        help="Path to fasta file with transcripts derived from official annotation")
    parser.add_argument('--annotation_transcripts', default=None, action=FullPaths,
                        help="Path to fasta file with full transcripts set")
    parser.add_argument('--annotation_proteins', default=None, action=FullPaths,
                        help="Path to fasta file with full proteins set")
    parser.add_argument('--external_gff', default=None, action=FullPaths,
                        help="Path to an external gff with gene models to be used")
    parser.add_argument('--busco_set', default=None, action=FullPaths,
                        help="Path to BUSCOs data set to be used for QA")
    parser.add_argument('--blast_qa_db', default=None, action=FullPaths,
                        help="Path to to a Blast proteins DB against which annotated proteins will be compared for QA")
    parser.add_argument('--busco_exe', action=FullPaths, required=True, help="Path to BUSCO executabl")
    parser.add_argument('--augustus_bin', action=FullPaths, required=True, help="Path to Augustus bin dir")
    parser.add_argument('--augustus_conf', action=FullPaths, required=True, help="Path to Augustus config file")
    parser.add_argument('--qa_report_script', action=FullPaths, required=True, help="Path to annotation QA script")
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
    genome_annotation_pipeline(**vars(args))
