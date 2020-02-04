"""
This pipeline constructs a pan genome
in the "map-to-pan" approach. It consists
of the following general steps:
1. Download fastq files from ena
2. Assemble reads from each sample into
   contigs
3. Iteratively map contigs to reference
   and adding novel sequences to create the
   non-reference section of the pan genome
4. Annotate the non-reference contigs
5. Align reads from each sample to reference
   + non-reference genes to determine gene
   presence/absence
6. Summarize and create PAV matrix
"""

import os
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir)
import sys
sys.path.append(utils_dir)
from snakemakeUtils import *
from collections import OrderedDict

def init():
    #load_info_file
    config['samples_info'] = SampleInfoReader.sample_table_reader(filename=config['samples_info_file'],
                delimiter='\t', key_name='sample', col_names=['ena_ref'])
init()
config['samples_info'] = OrderedDict(config['samples_info'])
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

localrules: all, prep_liftover_chunks_tsv, prep_annotation_chunks_tsv, prep_liftover_yaml, prep_annotation_yaml 

n_samples = len(config['samples_info'])
last_sample = list(config['samples_info'].keys())[-1]
last_sample_ena = config['samples_info'][last_sample]['ena_ref']
rule all:
    input:
        config["out_dir"] + "/all_samples/non_ref/non_redun_non_ref_contigs.fasta"

def get_sample(wildcards):
    return config['samples_info'][wildcards.sample]['ena_ref']


rule download_fastq:
    """
    Download reads data from ENA
    """
    output:
        config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_1.fastq.gz",
        config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_2.fastq.gz"
    params:
        sample_out_dir=config["out_dir"] + "/per_sample/{sample}/data",
        ena_ref=get_sample,
        download_script=utils_dir + '/ena-fast-download.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    #conda:
    #    CONDA_ENV_DIR + '/ena_download.yml'
    shell:
        """
        module load curl
        python {params.download_script} {params.ena_ref} --output_directory {params.sample_out_dir}
        """

rule quality_trimming:
    """
    Trim/remove low quality reads
    """
    input:
        r1=config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_1.fastq.gz",
        r2=config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_2.fastq.gz"
    output:
        r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_paired.fastq.gz",
        r1_unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_unpaired.fastq.gz",
        r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz",
        r2_unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_unpaired.fastq.gz"
    params:
        trimming_modules=config['trimming_modules'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/trimmomatic.yml'
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} {params.trimming_modules} -threads {params.ppn}
        """

rule merge_reads:
    """
    Merge read pairs to create long fragments
    """
    input:
        r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_paired.fastq.gz",
        r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz"
    output:
        config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.extendedFrags.fastq.gz",
        config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_1.fastq.gz",
        config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_2.fastq.gz"
    params:
        merge_out_dir=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}",
        merge_min_overlap=config['merge_min_overlap'],
        merge_max_mismatch_ratio=config['merge_max_mismatch_ratio'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/flash.yml'
    shell:
        """
        flash {input.r1_paired} {input.r2_paired} -d {params.merge_out_dir} -m {params.merge_min_overlap} -x {params.merge_max_mismatch_ratio} -z -t {params.ppn} -o {wildcards.ena_ref}
        """

rule combine_unpaired:
    """
    Combine unpaired R1 and R2 into one file
    """
    input:
        r1_unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_unpaired.fastq.gz",
        r2_unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_unpaired.fastq.gz"
    output:
        config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_clean_unpaired.fastq.gz"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        cat {input.r1_unpaired} {input.r2_unpaired} > {output}
        """

rule genome_assembly:
    """
    De novo assembly of reads into contigs
    """
    input:
        r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_1.fastq.gz",
        r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_2.fastq.gz",
        unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_clean_unpaired.fastq.gz",
        merged=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.extendedFrags.fastq.gz"
    output:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs.fasta",
    params:
        out_dir=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}",
        ppn=config['ppn'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/spades.yml'
    shell:
        """
        spades.py -o {params.out_dir} --pe1-1 {input.r1_paired} --pe1-2 {input.r2_paired} --pe1-m {input.merged} --pe1-s {input.unpaired} --threads {params.ppn}
        """

rule filter_contigs:
    """
    Discard contigs shorter than L or with coverage lower than C (L, C given by user)
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta"
    params:
        filter_script=utils_dir + '/filter_contigs.py',
        min_length=config['min_length'],
        min_coverage=config['min_coverage'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        python {params.filter_script} {input} {params.min_length} {params.min_coverage} {output}
        """

rule assembly_busco:
    """
    Run BUSCO on filtered contigs
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/run_BUSCO/short_summary_BUSCO.txt"
    params:
        assembly_dir=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}",
        busco_set=config['busco_set'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/busco.yml'
    shell:
        """
        cd {params.assembly_dir}
        run_busco -i {input} -o BUSCO -m genome -l {params.busco_set} -c {params.ppn} -f
        """

rule assembly_quast:
    """
    Run QUAST on filtered assembly to get assembly stats and QA
    """
    input:
        contigs=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta",
        r1=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_paired.fastq.gz",
        r2=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz"
    output:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/QUAST/report.html"
    params:
        out_dir=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/QUAST",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/quast.yml'
    shell:
        """
        quast {input.contigs} -o {params.out_dir} -t {params.ppn} -1 {input.r1} -2 {input.r2}
        """

#def recurse_sample(wcs):
#    n = int(wcs.n)
#    if n == 1:
#        ret = config['reference_genome']
#    elif n > 1:
#        ret = config["out_dir"] + '/all_samples/non_reference/%s_%s_%s/non_ref_contigs.fasta' % (wcs.sample, wcs.ena_ref, str(n-1))
#    else:
#        raise ValueError("loop numbers must be 1 or greater: received %s" % wcs.n)
#    return ret
#
#rule map_to_pan:
#    """
#    Collect the non-reference part of the
#    pan genome by iteratively mapping asembled
#    contigs to the reference
#    """
#    input:
#        ref=recurise_sample,
#        contigs=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta"
#    output:
#        map_sam=config["out_dir"] + '/all_samples/non_reference/{sample}_{ena_ref}_{n}/contigs_map.sam',
#        non_ref=config["out_dir"] + '/all_samples/non_reference/{sample}_{ena_ref}_{n}/non_ref_contigs.fasta'
#    params:
#        min_length=config['min_length'],
#        extract_insert_script=pipeline_dir + '/extract_insertions_from_sam.py',
#        queue=config['queue'],
#        priority=config['priority'],
#        logs_dir=LOGS_DIR,
#        ppn=config['ppn']
#    conda:
#        CONDA_ENV_DIR + '/map_to_pan.yml'
#    shell:
#        """
#        minimap2 -ax asm5 -t {params.ppn} {input.contigs} {input.ref} > {output.map_sam}
#        samtools view -f 4 tmp/map5.sam | awk 'length($10) > {params.min_length}' | samtools fasta - > {output.non_ref}
#        python {params.extract_insert_script} {output.map_sam} {params.min_length} >> {output.non_ref}
#        """

rule map_to_ref:
    """
    Map contigs from each genome assembly
    to reference in order to detect novel
    sequences
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/non_ref_{ena_ref}/contigs_filter_vs_ref.sam"
    params:
        ref_genome=config['reference_genome'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/minimap2.yml'
    shell:
        """
        minimap2 -ax asm5 -t {params.ppn} {params.ref_genome} {input} > {output}
        """

rule extract_unmapped:
    """
    Extract contigs not mapped to ref
    and convert to fasta
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/non_ref_{ena_ref}/contigs_filter_vs_ref.sam"
    output:
        config["out_dir"] + "/per_sample/{sample}/non_ref_{ena_ref}/contigs_filter_unmapped.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools fasta -f 4 {input} > {output}
        """

rule extract_inserts:
    """
    Extract large insert sequences
    from mapped contigs 
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/non_ref_{ena_ref}/contigs_filter_vs_ref.sam"
    output:
        config["out_dir"] + "/per_sample/{sample}/non_ref_{ena_ref}/contigs_filter_inserts.fasta"
    params:
        min_length=config['min_length'],
        extract_insert_script=pipeline_dir + '/extract_insertions_from_sam.py',        
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/pysam.yml'
    shell:
        """
        python {params.extract_insert_script} {input} {params.min_length} > {output}
        """

rule prep_non_ref:
    """
    Combine unmapped and insert sequences
    and add the sample name to each fasta
    record.
    """
    input:
        unmapped=config["out_dir"] + "/per_sample/{sample}/non_ref_{ena_ref}/contigs_filter_unmapped.fasta",
        insert=config["out_dir"] + "/per_sample/{sample}/non_ref_{ena_ref}/contigs_filter_inserts.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/non_ref_{ena_ref}/non_ref.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        cat {input.unmapped} {input.insert} | sed 's/>/>{wildcards.sample}_/' > {output}
        """

rule collect_all_non_ref:
    """
    Concat all (redundant) non-ref
    sequences into one fasta file
    """
    input:
        expand(config["out_dir"] + "/per_sample/{sample}/non_ref_{ena_ref}/non_ref.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()])
    output:
        config["out_dir"] + "/all_samples/non_ref/redun_non_ref_contigs.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        cat {input} > {output}
        """

rule remove_redundant:
    """
    Remove redundant non-ref sequences
    by clustering with CD-HIT and taking
    longest sequence from each cluster.
    """
    input:
        config["out_dir"] + "/all_samples/non_ref/redun_non_ref_contigs.fasta"
    output:
        config["out_dir"] + "/all_samples/non_ref/non_redun_non_ref_contigs.fasta"
    params:
        similarity_threshold=config['similarity_threshold'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn'],
    conda:
        CONDA_ENV_DIR + '/cd-hit.yml'
    shell:
        """
        cd-hit -i {input} -o {output} -c {params.similarity_threshold} -n 5 -M 0 -d 0 -T {params.ppn}
        """
