import os
import pandas as pd
import numpy as np

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

genomes_df = pd.read_csv(config['genomes_tsv'], sep='\t')
all_genomes = list(genomes_df['genome_name'])
if len(all_genomes) != len(set(all_genomes)):
    exit('Genome names must be unique!')
if not all(genomes_df['genome_type'].isin(['REF','HQ','LQ'])):
    exit('Invalid genome types found. Valid types: REF, HQ, LQ')
ref_genome = genomes_df.query('genome_type == "REF"')['genome_name'].values[0]
LQ_genomes = list(genomes_df.query('genome_type == "LQ"')['genome_name'].values)


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

def get_genome(wildcards):
    genome = wildcards.genome
    query = 'genome_name == @genome'
    fasta = genomes_df.query(query)['fasta'].values[0]
    return fasta

def get_mem_gb(wildcards, attempt):
    return attempt * 8
    
include: "rules/KMC.smk"
include: "rules/create_kmer_matrix.smk"
include: "rules/create_kinship_matrix.smk"
include: "rules/prepare_phenotypes.smk"
include: "rules/kmerGWAS.smk"
include: "rules/extract_significant_kmers.smk"
include: "rules/map_significant_kmers.smk"
include: "rules/filter_kmer_mapping.smk"
include: "rules/index_genome.smk"
include: "rules/extract_significant_contigs.smk"
include: "rules/map_significant_contigs_to_ref.smk"
include: "rules/filter_contig_mapping.smk"
include: "rules/project_kmers_to_ref.smk"
include: 'rules/create_final_results_table.smk'

wildcard_constraints:
    genome = '[A-Za-z0-9_\-]+',
    phenotype = '[A-Za-z0-9_\-]+'

rule all:
    input:
        expand(os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per_by_acc.tsv'), phenotype=phenotypes_list),
        expand(os.path.join(out_dir, 'all_samples', '{phenotype}', '{phenotype}_histogram.html'), phenotype=phenotypes_list),
        expand(os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per_mapping.tsv'), phenotype=phenotypes_list)

for r in workflow.rules:
    r.set_params(queue=config["queue"])
    r.set_params(priority=config["priority"])
