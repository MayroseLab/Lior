"""
This pipeline constructs a pan genome
in the "de novo" approach. It consists
of the following general steps:
1. Download fastq files from ena
2. Assemble reads from each sample into
   contigs
3. Annotate each assembly
4. Perform QA and filtration on annotations
5. Cluster genes from all samples to detect
   orthologs
6. Summarize and create PAV matrix
"""

import os
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir)
import sys
print(sys.version)
sys.path.append(utils_dir)
from snakemakeUtils import *

def init():
    #load_info_file
    config['samples_info'] = SampleInfoReader.sample_table_reader(filename=config['samples_info_file'],
                delimiter='\t', key_name='sample', col_names=['ena_ref'])
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

localrules: all, prep_liftover_chunks_tsv, prep_annotation_chunks_tsv, prep_liftover_yaml, prep_annotation_yaml 

rule all:
    input:
        expand(config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/run_BUSCO/short_summary_BUSCO.txt", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        expand(config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/QUAST/report.html", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        expand(config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/run_BUSCO/short_summary_BUSCO.txt", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()])

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

rule prep_chunks:
    """
    Divide filtered contigs into chunks for efficient parallel analysis
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta",
    output:
        config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/done"
    params:
        split_script=utils_dir + '/split_fasta_records.py',
        chunks=config['max_jobs'] - 1,
        out_dir=config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        python --version
        python {params.split_script} {input} {params.out_dir} {params.chunks}
        touch {output}
        """ 

rule prep_liftover_chunks_tsv:
    """
    Prepare TSV config for liftover run
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/done"
    output:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/chunks.tsv"
    params:
        chunks_dir=config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        echo "chunk\tpath" > {output}
        realpath {params.chunks_dir}/*.fasta | awk '{{n=split($0,a,"/"); sub(".yml","",a[n]); print a[n]"\t"$0}}' >> {output}
        """

rule prep_liftover_yaml:
    """
    Prepare yml config for liftover run
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/chunks.tsv"
    output:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/liftover.yml"
    params:
        liftover_dir=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}",
        templates_dir=config["liftover_config_templates"],
        liftover_transcripts=config['liftover_transcripts'],
        repeats_library=config['repeats_library'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        echo "name: MAKER_wrapper" >> {output}
        echo "chunks_info_file: {input}" >> {output}
        echo "out_dir: {params.liftover_dir}" >> {output}
        echo "config_templates: {params.templates_dir}" >> {output}
        echo "queue: {params.queue}" >> {output}
        echo "priority: {params.priority}" >> {output}
        echo "sample: {wildcards.sample}" >> {output}
        echo "logs_dir: {params.logs_dir}" >> {output}
        echo config_kv_pairs: est={params.liftover_transcripts} rmlib={params.repeats_library} >> {output}
        """

rule maker_liftover:
    """
    Run MAKER while only allowing liftover of reference transcripts
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/liftover.yml"
    output:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/maker.genes.gff"
    params:
        run_maker_in_chunks_snakefile=annotation_pipeline_dir + '/run_MAKER_in_chunks.snakefile',
        queue=config['queue'],
        jobs=config['max_jobs']//len(config['samples_info']),
        liftover_dir=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}",
        qsub_wrapper_script=utils_dir + '/pbs_qsub_snakemake_wrapper.py',
        priority=config['priority'],
        jobscript=utils_dir + '/jobscript.sh',
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        cd {params.liftover_dir}
        snakemake -s {params.run_maker_in_chunks_snakefile} --configfile {input} --cluster "python {params.qsub_wrapper_script}" -j {params.jobs} --latency-wait 60 --restart-times 3 --jobscript {params.jobscript}
            """

rule prep_annotation_chunks_tsv:
    """
    Prepare TSV config for annotation run
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/chunks.tsv"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/chunks.tsv"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cp {input} {output}
        """

rule prep_annotation_yaml:
    """
    Prepare yml config for annotation run
    """
    input:
        chunks_tsv=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/chunks.tsv",
        liftover_gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/maker.genes.gff"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/annotation.yml"
    params:
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}",
        templates_dir=config["annotation_config_templates"],
        liftover_transcripts=config['liftover_transcripts'],
        additional_transcripts=config['additional_transcripts'],
        proteins=config['proteins'],
        repeats_library=config['repeats_library'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        echo "name: MAKER_wrapper" >> {output}
        echo "chunks_info_file: {input.chunks_tsv}" >> {output}
        echo "out_dir: {params.annotation_dir}" >> {output}
        echo "config_templates: {params.templates_dir}" >> {output}
        echo "queue: {params.queue}" >> {output}
        echo "priority: {params.priority}" >> {output}
        echo "sample: {wildcards.sample}" >> {output}
        echo "logs_dir: {params.logs_dir}" >> {output}
        echo config_kv_pairs: est={params.liftover_transcripts},{params.additional_transcripts} protein={params.proteins} rmlib={params.repeats_library} pred_gff={input.liftover_gff} >> {output}
        """

rule maker_annotation:
    """
    Run MAKER using liftover results + evidence
    to obtain a full genome annotation
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/annotation.yml"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.fasta"
    params:
        run_maker_in_chunks_snakefile=annotation_pipeline_dir + '/run_MAKER_in_chunks.snakefile',
        queue=config['queue'],
        jobs=config['max_jobs']//len(config['samples_info']),
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}",
        qsub_wrapper_script=utils_dir + '/pbs_qsub_snakemake_wrapper.py',
        priority=config['priority'],
        jobscript=utils_dir + '/jobscript.sh',
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        cd {params.annotation_dir}
        snakemake -s {params.run_maker_in_chunks_snakefile} --configfile {input} --cluster "python {params.qsub_wrapper_script}" -j {params.jobs} --latency-wait 60 --restart-times 3 --jobscript {params.jobscript}
            """
rule filter_annotation:
    """
    Remove unreliable annotations
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter.fasta"
    params:
        filter_script=utils_dir + '/filter_by_aed.py',
        max_aed=config['max_aed'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        python {params.filter_script} {input} {params.max_aed} {output}
        """

rule annotation_busco:
    """
    Run BUSCO on filtered annotation proteins
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/run_BUSCO/short_summary_BUSCO.txt"
    params:
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}",
        busco_set=config['busco_set'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/busco.yml'
    shell:
        """
        cd {params.annotation_dir}
        run_busco -i {input} -o BUSCO -m proteins -l {params.busco_set} -c {params.ppn} -f
        """
