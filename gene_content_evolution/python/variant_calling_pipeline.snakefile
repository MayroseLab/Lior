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

def init():
    #load_info_file
    config['samples_info'] = SampleInfoReader.sample_table_reader(filename=config['samples_info_file'],
                delimiter='\t', key_name='sample', col_names=['ena_ref'])
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
        config["out_dir"] + "/all_variants_SNP.vcf"

def get_sample(wildcards):
    return config['samples_info'][wildcards.sample]['ena_ref']

def replace_extension(path,new_ext):
    return '.'.join([os.path.splitext(path)[0]] + [new_ext])

pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))

rule download_fastq:
    output:
        config["out_dir"] + "/per_sample/{sample}/{ena_ref}_1.fastq.gz",
        config["out_dir"] + "/per_sample/{sample}/{ena_ref}_2.fastq.gz"
    params:
        sample_out_dir=config["out_dir"] + "/per_sample/{sample}",
        ena_ref=get_sample,
        download_script=pipeline_dir + '/ena-fast-download.py',
        nodes=1,
        ppn=1,
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/curl.yml"
    shell:
        """
        python {params.download_script} {params.ena_ref} --output_directory {params.sample_out_dir}
        """


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
        R1=config["out_dir"] + "/per_sample/{sample}/{ENA_REF}_1.fastq.gz",
        R2=config["out_dir"] + "/per_sample/{sample}/{ENA_REF}_2.fastq.gz",
        index=config["reference_genome"] + '.bwt'
    output:
        config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.sam"
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
        config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.sam"
    output:
        config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.sort.bam"
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
        config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.sort.bam"
    output:
        md_bam=config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.sort.md.bam",
        md_stats=config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.sort.md.stats"
    params:
        nodes=1,
        ppn=1,
        ram=config['gatk_ram'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/gatk.yml"
    shell:
        """
        gatk --java-options "-Xmx{params.ram}" MarkDuplicates -I {input} -O {output.md_bam} -M {output.md_stats}
        """

rule index_bam:
    input:
        config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.sort.md.bam"
    output:
        config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.sort.md.bam.bai"
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
        replace_extension(config["reference_genome"],'dict')
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
        md_bam=config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.sort.md.bam",
        ref_dict=replace_extension(config["reference_genome"],'dict'),
        ref_faidx=config["reference_genome"] + '.fai',
        bai=config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.sort.md.bam.bai"
    output:
        config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.g.vcf"
    params:
        reference_genome=config["reference_genome"],
        nodes=1,
        ppn=1,
        ram=config['gatk_ram'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        "variant_calling_pipeline_conda_envs/gatk.yml"
    shell:
        """
        gatk --java-options "-Xmx{params.ram}" HaplotypeCaller -I {input.md_bam} -O {output} -R {params.reference_genome} -ERC GVCF
        """

rule combine_GVCFs:
    input:
        expand(config["out_dir"] + "/per_sample/{sample}/{sample}_{ENA_REF}_vs_ref.g.vcf", zip, sample=config['samples_info'].keys(), ENA_REF=[x['ena_ref'] for x in config['samples_info'].values()])
    output:
        config["out_dir"] + "/all_variants.g.vcf"
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
        files=`echo {input} | sed 's/ / -V /g'`
        gatk CombineGVCFs -V $files -O {output} -R {params.reference_genome} 
        """

rule genotype_samples:
    input:
        config["out_dir"] + "/all_variants.g.vcf"
    output:
        config["out_dir"] + "/all_variants.vcf"
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
        gatk GenotypeGVCFs -V {input} -O {output} -R {params.reference_genome}
        """

rule extract_SNPs:
    input:
        config["out_dir"] + "/all_variants.vcf"
    output:
        config["out_dir"] + "/all_variants_SNP.vcf"
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
        gatk SelectVariants --select-type-to-include SNP -V {input} -O {output}
        """
