from Bio import SeqIO
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir)
import sys
sys.path.append(utils_dir)
from snakemakeUtils import *

def init():
    #load_info_file
    config['pan_genomes'] = SampleInfoReader.sample_table_reader(filename=config['pan_genomes_info'],
                delimiter='\t', key_name='pan_genome_name', col_names=['proteins_fasta','pav_tsv'])
    assert len(config['pan_genomes']) == 2, "Exactly two pan genomes should be provided"

init()
LOGS_DIR = config['out_dir'] + "/logs"
CONDA_ENV_DIR = pipeline_dir + "/conda_env"
annotation_pipeline_dir = os.path.dirname(pipeline_dir) + '/annotation_pipeline'

onstart:
    write_config_file(config)
    if not os.path.isdir(LOGS_DIR):
        os.mkdir(LOGS_DIR)

onsuccess:
    print("%s pipeline finished, no error" % config['name'])
    shell("cat {log} >> %s/run_log.txt" % config["out_dir"])

onerror:
    print("%s pipeline failed" % config['name'])
    shell("cat {log} >> %s/run_log.txt" % config["out_dir"])

#------------------------------------
#                RULES              |
#------------------------------------

rule all:
    input:
        os.path.join(config["out_dir"], 'RBBH.tsv')

def get_pan_genome(wildcards):
    return config['pan_genomes'][wildcards.PG]['proteins_fasta']

rule extract_non_ref_PG:
    """
    Extract non-ref proteins into a new fasta
    """
    input:
        get_pan_genome
    output:
        os.path.join(config["out_dir"], "{PG}_nonref.fasta")
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    run:
        nonref_records = []
        for rec in SeqIO.parse(input[0], 'fasta'):
            if rec.id.startswith('PanGene'):
                nonref_records.append(rec)
        SeqIO.write(nonref_records, output[0], 'fasta')

rule make_blast_db:
    """
    Create blast DB for non-ref proteins of each pan genome
    """
    input:
        os.path.join(config["out_dir"], "{PG}_nonref.fasta")
    output:
        os.path.join(config["out_dir"], "{PG}_nonref.fasta.phr")
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/blast.yml'
    shell:
        """
        makeblastdb -dbtype prot -in {input} -input_type fasta
        """

rule blast_non_ref:
    """
    blast non ref proteins from one PG against the other
    """
    input:
        pg1_fasta=os.path.join(config["out_dir"], "{PG1}_nonref.fasta"),
        pg2_fasta=os.path.join(config["out_dir"], "{PG2}_nonref.fasta"),
        pg2_db=os.path.join(config["out_dir"], "{PG2}_nonref.fasta.phr")
    output:
        os.path.join(config["out_dir"], "{PG1}_nonref_vs_{PG2}_nonref.blast6")
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/blast.yml'
    shell:
        """
        blastp -query {input.pg1_fasta} -db {input.pg2_fasta} -out {output} -max_target_seqs 1 -outfmt 6
        """

pg_names = sorted(list(config['pan_genomes'].keys()))
pg1 = pg_names[0]
pg2 = pg_names[1]

rule find_RBBH:
    """
    Find reciprocal best blast hits between
    non-ref proteins of the two pan genomes
    """
    input:
        fw=os.path.join(config["out_dir"], pg1 + "_nonref_vs_" + pg2 + "_nonref.blast6"),
        rev=os.path.join(config["out_dir"], pg2 + "_nonref_vs_" + pg1 + "_nonref.blast6")
        #expand(os.path.join(config["out_dir"], "{PG1}_nonref_vs_{PG2}_nonref.blast6"), zip, PG1=sorted(list(config['pan_genomes'].keys())), PG2=sorted(list(config['pan_genomes'].keys())))
    output:
        os.path.join(config["out_dir"], 'RBBH.tsv')
    params:
        find_RBBH_script=os.path.join(pipeline_dir, 'find_RBBH.py'), 
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        python {params.find_RBBH_script} {input.fw} {input.rev} {output}
        """




