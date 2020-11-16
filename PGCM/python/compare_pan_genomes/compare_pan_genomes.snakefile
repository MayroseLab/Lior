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
    assert len(config['pan_genomes']) == 3, "Exactly three pan genomes should be provided"
    assert "TRUE_PG" in config['pan_genomes'], "True pan-genome must be included and named TRUE_PG"

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

localrules: all

rule all:
    input:
        os.path.join(config["out_dir"], 'report.html')

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
        blastp -query {input.pg1_fasta} -db {input.pg2_fasta} -out {output} -max_target_seqs 5 -outfmt 6
        """

pg_names = sorted(list(config['pan_genomes'].keys()))
pg1 = pg_names[0]
pg2 = pg_names[1]
true_pg = 'TRUE_PG'

rule find_matches:
    """
    Find the maximum weight matching between pan genome proteins
    """
    input:
        pg1_fasta=os.path.join(config["out_dir"], "{PG1}_nonref.fasta"),
        pg2_fasta=os.path.join(config["out_dir"], "{PG2}_nonref.fasta"),
        fw=os.path.join(config["out_dir"], "{PG1}_nonref_vs_{PG2}_nonref.blast6"),
        rev=os.path.join(config["out_dir"], "{PG2}_nonref_vs_{PG1}_nonref.blast6"),
    output:
        os.path.join(config["out_dir"], '{PG1}_vs_{PG2}_max_weight_matches.tsv')
    params:
        find_matches_script=os.path.join(pipeline_dir, 'match_non_ref.py'),
        pg1_name=lambda w: w.PG1,
        pg2_name=lambda w: w.PG2,
        min_bitscore=config['min_bitscore'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/match_non_ref.yml'
    shell:
        """
        python {params.find_matches_script} {input.pg1_fasta} {input.pg2_fasta} {input.fw} {input.rev} {output} --min_weight {params.min_bitscore} --set1_name {params.pg1_name} --set2_name {params.pg2_name}
        """

rule create_report_nb:
    """
    Create jupyter notebook of comparison report
    """
    input:
        pg1_pav=config['pan_genomes'][pg1]['pav_tsv'],
        pg2_pav=config['pan_genomes'][pg2]['pav_tsv'],
        true_pg_pav=config['pan_genomes']['TRUE_PG']['pav_tsv'],
        pg1_vs_pg2_matches=os.path.join(config["out_dir"], '{}_vs_{}_max_weight_matches.tsv'.format(pg1,pg2,true_pg)),
        pg1_vs_true_matches=os.path.join(config["out_dir"], '{}_vs_{}_max_weight_matches.tsv'.format(pg1,true_pg,true_pg)),
        pg2_vs_true_matches=os.path.join(config["out_dir"], '{}_vs_{}_max_weight_matches.tsv'.format(pg2,true_pg,true_pg))
    output:
        os.path.join(config["out_dir"], 'report.ipynb')
    params:
        nb_template=os.path.join(pipeline_dir, 'report_template.ipynb'),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        sed -e 's@<PG1_PAV>@{input.pg1_pav}@' -e 's@<PG2_PAV>@{input.pg2_pav}@' -e 's@<TRUE_PAV>@{input.true_pg_pav}@' -e 's@<PG1_VS_PG2_NON_REF_MATCHES>@{input.pg1_vs_pg2_matches}@' -e 's@<PG1_VS_TRUE_NON_REF_MATCHES>@{input.pg1_vs_true_matches}@' -e 's@<PG2_VS_TRUE_NON_REF_MATCHES>@{input.pg2_vs_true_matches}@' -e 's@<PG1_NAME>@%s@' -e 's@<PG2_NAME>@%s@' {params.nb_template} > {output} 
        """ %(pg1, pg2)

rule create_report_html:
    """
    Convert notebook to HTML report
    """
    input:
        os.path.join(config["out_dir"], 'report.ipynb')
    output:
        os.path.join(config["out_dir"], 'report.html')
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/jupyter.yml'
    shell:
        """
        jupyter nbconvert {input} --output {output} --no-prompt --no-input --execute --NotebookClient.timeout=-1 --ExecutePreprocessor.timeout=-1 --to html
        """
