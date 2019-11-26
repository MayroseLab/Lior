"""
This pipeline performs end-to-end variant
calling from sequencing reads, including
the following steps:
1. Download fastq data
2. Reads cleanup and pre-processing
3. Reads alignment
4. Variant calling
5. Variant filtration

It can run on multiple samples and will
produce a VCF file per sample.
See the variant_calling_pipeline.conf.yml
file for input parameters
"""

from snakemakeUtils import *
import os
from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
NCBI = NCBIRemoteProvider(email="liorglic@mail.tau.ac.il")

def init():
    #load_info_file
    config['samples_info'] = SampleInfoReader.sample_table_reader(filename=config['samples_info_file'],
                delimiter='\t', key_name='sample', col_names=['R1','R2'])
init()

LOGS_DIR = config['out_dir'] + "/logs"

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
        config["out_dir"] + "/all_variants.vcf"

def get_sample_R1(wildcards):
    return config['samples_info'][wildcards.sample]['r1']
def get_sample_R2(wildcards):
    return config['samples_info'][wildcards.sample]['r2']


#rule download_fastq:
#    input:
#
#    output:
#        config["out_dir"] + "/per_sample/{sample}/{ENA_REF}_1.fastq.gz",
#        config["out_dir"] + "/per_sample/{sample}/{ENA_REF}_2.fastq.gz"
#    params:
#        sample_out_dir=config["out_dir"] + "/per_sample/{sample}",
#        nodes=1,
#        ppn=1,
#        queue=config['queue'],
#        priority=config['priority'],
#        logs_dir=LOGS_DIR
#    conda:
#        "variant_calling_pipeline_conda_envs/curl.yml"
#    shell:
#        """
#        python ena-fast-download.py {params.ENA_REF} --output_directory {params.sample_out_dir}
#        """


rule index_ref_genome:
    input:
        config["reference_genome"]
    output:
        config["reference_genome"] + '.bwt'
    params:
        nodes=1,
        ppn=1,
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/bwa.yml"
    shell:
        """
        bwa index {input}
        """

rule map_to_ref:
    input:
        R1=get_sample_R1,
        R2=get_sample_R2,
        index=config["reference_genome"] + '.bwt'
    output:
        config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.sam"
    params:
        reference_genome=config["reference_genome"],
        nodes=1,
        ppn=config['ppn'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/bwa.yml"
    shell:
        """
        bwa mem {params.reference_genome} {input.R1} {input.R2} -t {params.ppn} -o {output} -R '@RG\\tID:1\\tLB:READS\\tPL:ILLUMINA\\tSM:{wildcards.sample}\\tPU:1'
        """

rule sam_to_sorted_bam:
    input:
        config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.sam"
    output:
        config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.sort.bam"
    params:
        nodes=1,
        ppn=config['ppn'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/samtools.yml"
    shell:
        """
        samtools view -@ {params.ppn} -bh {input} | samtools sort -@ {params.ppn} - -o {output}
        """

rule mark_duplicates:
    input:
        config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.sort.bam"
    output:
        md_bam=config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.sort.md.bam",
        md_stats=config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.sort.md.stats"
    params:
        nodes=1,
        ppn=1,
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/gatk.yml"
    shell:
        """
        gatk MarkDuplicates -I {input} -O {output.md_bam} -M {output.md_stats}
        """

rule index_bam:
    input:
        config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.sort.md.bam"
    output:
        config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.sort.md.bam.bai"
    params:
        nodes=1,
        ppn=1,
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/samtools.yml"
    shell:
        """
        samtools index {input}
        """

rule create_ref_dict:
    input:
        config["reference_genome"]
    output:
        config["reference_genome"].replace('fasta','dict')
    params:
        nodes=1,
        ppn=1,
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/gatk.yml"
    shell:
        """
        gatk CreateSequenceDictionary -R {input} -O {output}
        """

rule create_ref_faidx:
    input:
        config["reference_genome"]
    output:
        config["reference_genome"] + '.fai'
    params:
        nodes=1,
        ppn=1,
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/samtools.yml"
    shell:
        """
        samtools faidx {input}
        """

rule call_variants:
    input:
        md_bam=config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.sort.md.bam",
        ref_dict=config["reference_genome"].replace('fasta','dict'),
        ref_faidx=config["reference_genome"] + '.fai',
        bai=config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.sort.md.bam.bai"
    output:
        config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.g.vcf"
    params:
        reference_genome=config["reference_genome"],
        nodes=1,
        ppn=1,
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/gatk.yml"
    shell:
        """
        gatk HaplotypeCaller -I {input.md_bam} -O {output} -R {params.reference_genome} -ERC GVCF
        """

rule combine_GVCFs:
    input:
        expand(config["out_dir"] + "/per_sample/{sample}/{sample}_vs_ref.g.vcf", sample=config['samples_info'])
    output:
        config["out_dir"] + "/all_variants.g.vcf"
    params:
        reference_genome=config["reference_genome"],
        nodes=1,
        ppn=config['ppn'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/gatk.yml"
    shell:
        """
        gatk CombineGVCFs -V {input} -O {output} -R {params.reference_genome} 
        """

rule genotype_samples:
    input:
        config["out_dir"] + "/all_variants.g.vcf"
    output:
        config["out_dir"] + "/all_variants.vcf"
    params:
        reference_genome=config["reference_genome"],
        nodes=1,
        ppn=config['ppn'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/gatk.yml"
    shell:
        """
        gatk GenotypeGVCFs -V {input} -O {output} -R {params.reference_genome}
        """
