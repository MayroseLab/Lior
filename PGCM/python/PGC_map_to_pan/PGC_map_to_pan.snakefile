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
         config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv"

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

rule prep_annotation_chunks:
    """
    Divide non-ref contigs into chunks for efficient parallel analysis
    """
    input:
        config["out_dir"] + "/all_samples/non_ref/non_redun_non_ref_contigs.fasta"
    output:
        config["out_dir"] + "/all_samples/non_ref/chunks/chunks.done"
    params:
        split_script=utils_dir + '/split_fasta_records.py',
        chunks=config['max_jobs'] - 1,
        out_dir=config["out_dir"] + "/all_samples/non_ref/chunks",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        python {params.split_script} {input} {params.out_dir} {params.chunks}
        touch {output}
        """ 

rule prep_annotation_chunks_tsv:
    """
    Prepare TSV config for liftover run
    """
    input:
        config["out_dir"] + "/all_samples/non_ref/chunks/chunks.done"
    output:
        config["out_dir"] + "/all_samples/non_ref/chunks/chunks.tsv"
    params:
        chunks_dir=config["out_dir"] + "/all_samples/non_ref/chunks",
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

rule prep_annotation_yaml:
    """
    Prepare yml config for annotation run
    """
    input:
        chunks_tsv=config["out_dir"] + "/all_samples/non_ref/chunks/chunks.tsv"
    output:
        config["out_dir"] + "/all_samples/annotation/annotation.yml"
    params:
        annotation_dir=config["out_dir"] + "/all_samples/annotation",
        templates_dir=config["annotation_config_templates"],
        transcripts=config['transcripts'],
        proteins=config['proteins'],
        repeats_library=config['repeats_library'],
        augustus_species=config['augustus_species'],
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
        echo "sample: non_ref_contigs" >> {output}
        echo "logs_dir: {params.logs_dir}" >> {output}
        echo config_kv_pairs: est={params.transcripts} protein={params.proteins} rmlib={params.repeats_library} augustus_species={params.augustus_species} >> {output}
        """

rule maker_annotation:
    """
    Run MAKER on non-ref contigs
    """
    input:
        config["out_dir"] + "/all_samples/annotation/annotation.yml"
    output:
        config["out_dir"] + "/all_samples/annotation/maker.proteins.fasta",
        config["out_dir"] + "/all_samples/annotation/maker.genes.gff"
    params:
        run_maker_in_chunks_snakefile=annotation_pipeline_dir + '/run_MAKER_in_chunks.snakefile',
        queue=config['queue'],
        jobs=config['max_jobs'],
        annotation_dir=config["out_dir"] + "/all_samples/annotation",
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

rule rename_genes:
    """
    Assign genes short, unique names (gff and fasta).
    Names consist of the genome name and a unique ID.
    """
    input:
        fasta=config["out_dir"] + "/all_samples/annotation/maker.proteins.fasta",
        gff=config["out_dir"] + "/all_samples/annotation/maker.genes.gff"
    output:
        config["out_dir"] + "/all_samples/annotation/maker.proteins.raw_names.fasta",
        config["out_dir"] + "/all_samples/annotation/maker.genes.raw_names.gff",
        config["out_dir"] + "/all_samples/annotation/gff.map",
        config["out_dir"] + "/all_samples/annotation/fasta.map"
    params:
        out_dir=config["out_dir"] + "/all_samples/annotation/",
        create_fasta_map_script=os.path.join(utils_dir,"create_fasta_map.py"),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        maker_map_ids --prefix PanGene --justify 1 --iterate 1 {input.gff} > {params.out_dir}/gff.map
        cp {input.gff} {params.out_dir}/maker.genes.raw_names.gff
        map_gff_ids {params.out_dir}/gff.map {input.gff}
        python {params.create_fasta_map_script} {input.gff} > {params.out_dir}/fasta.map
        cp {input.fasta} {params.out_dir}/maker.proteins.raw_names.fasta
        map_fasta_ids {params.out_dir}/fasta.map {input.fasta}
        """

rule filter_annotation:
    """
    Remove unreliable annotations
    """
    input:
        fasta=config["out_dir"] + "/all_samples/annotation/maker.proteins.fasta",
        fasta_map=config["out_dir"] + "/all_samples/annotation/fasta.map"
    output:
        config["out_dir"] + "/all_samples/annotation/maker.proteins_filter.fasta"
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
        python {params.filter_script} {input.fasta} {params.max_aed} {output}
        """

rule prevent_duplicate_names:
    """
    If for any reason the filtered
    proteins contain duplicate names,
    prevent this by renaming. This
    helps prevent problems in next steps
    """
    input:
        config["out_dir"] + "/all_samples/annotation/maker.proteins_filter.fasta"
    output:
        config["out_dir"] + "/all_samples/annotation/maker.proteins_filter_nodupl.fasta"
    params:
        dupl_script=utils_dir + '/prevent_duplicate_names.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.dupl_script} {input} {output}
        """

rule remove_redundant_proteins:
    """
    Remove redundant non-ref proteins
    by clustering with CD-HIT and taking
    longest sequence from each cluster.
    """
    input:
        config["out_dir"] + "/all_samples/annotation/maker.proteins_filter_nodupl.fasta"
    output:
        config["out_dir"] + "/all_samples/annotation/non_redun_maker.proteins_filter_nodupl.fasta"
    params:
        similarity_threshold=config['similarity_threshold_proteins'],
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

rule match_gff:
    """
    Filter genes gff according to
    the genes remaining in proteins
    fasta after all filtrations.
    """
    input:
        prot=config["out_dir"] + "/all_samples/annotation/non_redun_maker.proteins_filter_nodupl.fasta",
        gff=config["out_dir"] + "/all_samples/annotation/maker.genes.gff",
    output:
        gff=config["out_dir"] + "/all_samples/annotation/non_redun_maker.genes_filter_nodupl.gff",
        filter_list=config["out_dir"] + "/all_samples/annotation/filter.list"
    params:
        filter_gff_script=utils_dir + '/filter_gff_by_id_list.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        grep '>' {input.prot} | sed 's/>\(.*\)-R.*/\\1/' > {output.filter_list}
        grep '>' {input.prot} | sed 's/>\(.*-R[1-9]*\).*/\\1/' >> {output.filter_list}
        python {params.filter_gff_script} {input.gff} {output.filter_list} {output.gff}
        """

rule remove_ref_alt_splicing:
    """
    In case the reference annotation
    contains genes with multiple mRNAs,
    only keep the longest transcript.
    """
    input:
        config['reference_annotation']
    output:
        config["out_dir"] + "/all_samples/ref/" + config['ref_name'] + '_longest_trans.gff'
    params:
        longest_trans_script=utils_dir + '/remove_alt_splicing_from_gff.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.longest_trans_script} {input} {output}
        """

rule create_pan_genome:
    """
    Create pan genome nucleotide
    and protein fasta files + gff
    as ref + non-ref
    """
    input:
        non_ref_contigs=config["out_dir"] + "/all_samples/non_ref/non_redun_non_ref_contigs.fasta",
        non_ref_proteins=config["out_dir"] + "/all_samples/annotation/non_redun_maker.proteins_filter_nodupl.fasta",
        non_ref_gff=config["out_dir"] + "/all_samples/annotation/non_redun_maker.genes_filter_nodupl.gff",
        ref_gff=config["out_dir"] + "/all_samples/ref/" + config['ref_name'] + '_longest_trans.gff'
    output:
        pan_genome=config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta",
        pan_proteome=config["out_dir"] + "/all_samples/pan_genome/pan_proteome.fasta",
        pan_genes=config["out_dir"] + "/all_samples/pan_genome/pan_genes.gff"
    params:
        ref_genome=config['reference_genome'],
        ref_proteins=config['reference_proteins'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        cat {params.ref_genome} {input.non_ref_contigs} > {output.pan_genome}
        cat {input.ref_gff} {input.non_ref_gff} > {output.pan_genes}
        cat {params.ref_proteins} {input.non_ref_proteins} > {output.pan_proteome}
        """

rule index_pan_genome:
    """
    Index pan genome for BWA runs
    """
    input:
        config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta"
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta.bwt"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/bwa.yml'
    shell:
        """
        bwa index {input}
        """

rule map_reads_to_pan:
    """
    Map reads from each sample to
    the pan genome.
    """
    input:
        r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_paired.fastq.gz",
        r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz",
        pan_genome=config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta",
        pan_genome_index=config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta.bwt"
    output:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sam"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/bwa.yml'
    shell:
        """
        bwa mem -t {params.ppn} {input.pan_genome} {input.r1_paired} {input.r2_paired} > {output}
        """

rule sam_to_sorted_bam:
    """
    Convert SAM output to sorted BAM
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sam"
    output:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sort.bam"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools view -@ {params.ppn} -bh {input} | samtools sort -@ {params.ppn} - -o {output}
        """

rule index_bam:
    input:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sort.bam"
    output:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sort.bam.bai"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools index {input}
        """

rule detect_gene_loss:
    """
    Run SGSGeneLoss on reads mapping result
    to detect gene losses per sample
    """
    input:
        bam=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sort.bam",
        bai=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sort.bam.bai",
        gff=config["out_dir"] + "/all_samples/pan_genome/pan_genes.gff"
    output:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/graph.csv"
    params:
        SGSGeneLossJar=config['SGSGeneLossJar'],
        bamDir=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/",
        bamList="{ena_ref}_map_to_pan.sort.bam",
        outDir=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/",
        minCov=config['minCov'],
        lostCutoff=config['lostCutoff'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        java -Xmx4g -jar {params.SGSGeneLossJar} bamPath={params.bamDir} bamFileList={params.bamList} gffFile={input.gff} outDirPath={params.outDir} chromosomeList=all minCov={params.minCov} lostCutoff={params.lostCutoff}
        """

rule unify_gene_loss_chromosomes:
    """
    Per sample, create one csv for all
    chromosomes by unifying all .excov
    files.
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/graph.csv"
    output:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{sample}.all.PAV"
    params:
        wdir=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        echo "chromosome,ID,is_lost,start_position,end_postion,frac_exons_covered,frac_gene_covered,ave_cov_depth_exons,cov_cat,ave_cove_depth_gene" > {output}
        grep -h -v is_lost {params.wdir}/*.excov >> {output}
        """

rule create_pan_PAV_matrix:
    """
    Create a unified matrix of gene PAV
    across all samples
    """
    input:
        expand(config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{sample}.all.PAV", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()])
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv"
    params:
        create_PAV_matrix_script=pipeline_dir + '/create_PAV_matrix.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/pandas.yml'
    shell:
        """
        python {params.create_PAV_matrix_script} {input} {output}
        """
