import os
import pandas as pd

out_dir = config['out_dir']
logs_dir = os.path.join(out_dir, 'logs')
envs_dir = os.path.join(workflow.basedir, 'envs')
scripts_dir = os.path.join(workflow.basedir, 'scripts')

samples_df = pd.read_csv(config['samples_tsv'], sep='\t')
samples = list(samples_df['sample'])
if len(samples) != len(set(samples)):
    exit('Sample names must be unique!')

phenotypes_df = pd.read_csv(config['phenotypes_tsv'], sep='\t', index_col=0)
assert set(samples_df['sample']) == set(phenotypes_df.index), "Different samples found in samples and phenotypes tables. Terminating."
phenotypes_list = list(phenotypes_df.columns)

def get_reads(wildcards):
    sample = wildcards.sample
    reads = samples_df.query('sample == @sample')['reads'].values[0]
    reads_file_list = reads.split(',')
    extensions = set()
    for f in reads_file_list:
        if not os.path.isfile(f):
            exit('File %s for sample %s does not exist' %(f, sample))
    return reads_file_list

def get_file_type(wildcards):
    sample = wildcards.sample
    file_type = samples_df.query('sample == @sample')['file_type'].values[0]
    if file_type not in ['fa', 'fq']:
        exit('Unknown file type %s for sample {wildcards.sample} - must be fa or fq')
    return file_type
    

include: "rules/KMC.smk"
include: "rules/create_kmer_matrix.smk"
include: "rules/create_kinship_matrix.smk"
include: "rules/prepare_phenotypes.smk"
include: "rules/kmerGWAS.smk"

rule all:
    input:
        expand(os.path.join(out_dir, 'all_samples', '{phenotype}', 'GWAS', 'kmers/output/phenotype_value.assoc.txt'), phenotype=phenotypes_list),
        expand(os.path.join(out_dir, 'all_samples', '{phenotype}', '{phenotype}_histogram.html'), phenotype=phenotypes_list)

for r in workflow.rules:
    r.set_params(queue=config["queue"])
    r.set_params(priority=config["priority"])
