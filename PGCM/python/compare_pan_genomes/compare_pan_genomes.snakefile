from Bio import SeqIO
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir)
import sys
sys.path.append(utils_dir)
from snakemakeUtils import *
import pandas as pd

#load_info_file
pan_genomes = pd.read_table(config['pan_genomes_info']).set_index("pan_genome_name", drop=False)
assert pan_genomes.shape[0] == 3, "Exactly three pan genomes should be provided"
assert "TRUE_PG" in pan_genomes.index, "True pan-genome must be included and named TRUE_PG"

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
        os.path.join(config["out_dir"], 'report.html'),
        os.path.join(config["out_dir"], 'discrepancies.tsv')

def get(wildcards, what):
    sample_dir = pan_genomes.loc[wildcards.PG, 'path']
    if what == 'prot':
        return sample_dir + '/all_samples/pan_genome/pan_proteome.fasta'
    elif what == 'pav':
        return sample_dir + '/all_samples/pan_genome/pan_PAV.tsv'

rule extract_non_ref_PG:
    """
    Extract non-ref proteins into a new fasta
    """
    input:
        lambda wc: get(wc, what='prot')
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

pg1 = pan_genomes.index[0]
pg2 = pan_genomes.index[1]
true_pg = pan_genomes.index[2]

rule create_report_nb:
    """
    Create jupyter notebook of comparison report
    """
    input:
        pg1_pav=pan_genomes.loc[pg1]['path'] + '/all_samples/pan_genome/pan_PAV.tsv',
        pg2_pav=pan_genomes.loc[pg2]['path'] + '/all_samples/pan_genome/pan_PAV.tsv',
        true_pg_pav=pan_genomes.loc[true_pg]['path'] + '/all_samples/pan_genome/pan_PAV.tsv',
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
        report=os.path.join(config["out_dir"], 'report.html'),
        discrep=os.path.join(config["out_dir"], 'discrepancies.tsv')
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/jupyter.yml'
    shell:
        """
        jupyter nbconvert {input} --output {output.report} --no-prompt --no-input --execute --NotebookClient.timeout=-1 --ExecutePreprocessor.timeout=-1 --to html
        """
