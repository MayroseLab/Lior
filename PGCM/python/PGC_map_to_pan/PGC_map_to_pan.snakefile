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
    # load HQ genomes info file
    config['hq_info'] = SampleInfoReader.sample_table_reader(filename=config['hq_genomes_info_file'],
                delimiter='\t', key_name='sample', col_names=['annotation_gff','genome_fasta', 'proteins_fasta'])
    # ensure not duplicate sample names exist
    all_names = list(config['samples_info'].keys()) + list(config['hq_info'].keys())
    assert len(all_names) == len(set(all_names)), "Can't use duplicate sample names!"

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

localrules: all, prep_liftover_chunks_tsv, prep_annotation_chunks_tsv, prep_liftover_yaml, prep_annotation_yaml, collect_HQ_mapping_bams, collect_LQ_mapping_bams

n_samples = len(config['samples_info'])
last_sample = list(config['samples_info'].keys())[-1]
last_sample_ena = config['samples_info'][last_sample]['ena_ref']
rule all:
    input:
        config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv",
        config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta",
        config["out_dir"] + "/all_samples/pan_genome/pan_proteome.fasta",
        config["out_dir"] + "/all_samples/pan_genome/pan_genes.gff"

def get_sample(wildcards):
    return config['samples_info'][wildcards.sample]['ena_ref']

def get_hq_sample_gff(wildcards):
    return config['hq_info'][wildcards.sample]['annotation_gff']

def get_hq_sample_genome(wildcards):
    return config['hq_info'][wildcards.sample]['genome_fasta']

def get_hq_sample_proteins(wildcards):
    return config['hq_info'][wildcards.sample]['proteins_fasta']

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

rule simplify_ref_gff_ID:
    """
    Edit ID (and Parent) attributes of ref gff
    To make them simpler (if needed)
    """
    input:
        config['reference_annotation']
    output:
        config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_simp.gff'
    params:
        simplify_gff_script=utils_dir + '/simplify_gff_id.py',
        simp_function=config['id_simplify_function'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.simplify_gff_script} {input} {output} "{params.simp_function}"
        """

rule remove_ref_alt_splicing:
    """
    In case the reference annotation
    contains genes with multiple mRNAs,
    only keep the longest transcript.
    """
    input:
        config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_simp.gff'
    output:
        gff=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans_simp.gff',
        table=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans_simp.gff.gene_to_mRNA'
    params:
        longest_trans_script=utils_dir + '/remove_alt_splicing_from_gff.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.longest_trans_script} {input} {output.gff}
        """

rule iterative_map_to_pan_HQ:
    """
    Iteratively create HQ pan genome
    by mapping HQ samples to ref and
    adding novel sequences and genes.
    """
    input:
       samples=config['hq_genomes_info_file'],
       ref_genome=config['reference_genome'],
       ref_gff=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans_simp.gff',
       ref_proteins=config['reference_proteins']
    output:
       config["out_dir"] + "/HQ_samples/HQ_pan/pan_genome.fasta",
       config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes.gff",
       config["out_dir"] + "/HQ_samples/HQ_pan/pan_proteins.fasta"
    params:
        map_to_pan_script=pipeline_dir + '/iterative_map_to_pan.py',
        min_len=config['min_length'],
        out_dir=config["out_dir"] + "/HQ_samples/HQ_pan",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/iterative_map_to_pan.yml'
    shell:
        """
        python {params.map_to_pan_script} {input.ref_genome} {input.ref_gff} {input.ref_proteins} {input.samples} {params.out_dir} --cpus {params.ppn} --min_len {params.min_len}
        """

rule remove_alt_splicing_from_HQ_pan:
    """
    Remove splice variants from HQ pan gff.
    Only keep longest variant. Then match
    the proteins fasta
    """
    input:
        gff=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes.gff",
        prot=config["out_dir"] + "/HQ_samples/HQ_pan/pan_proteins.fasta"
    output:
        gff=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes_no_alt_splicing.gff3",
        prot=config["out_dir"] + "/HQ_samples/HQ_pan/pan_proteins_no_alt_splicing.fasta",
        mrna_to_gene=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes_no_alt_splicing.gff3" + '.gene_to_mRNA'
    params:
        longest_trans_script=utils_dir + '/remove_alt_splicing_from_gff.py',
        filter_fasta_script=utils_dir + '/filter_fasta_by_gff.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.longest_trans_script} {input.gff} {output.gff}
        python {params.filter_fasta_script} {output.gff} {input.prot} {output.prot} mRNA ID
        """

rule prep_tsv_for_LQ_samples:
    """
    Prepare the TSV file required as
    input for creating the pan genome
    from assembled LQ samples
    """
    input:
        expand(config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/scaffolds.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()])
    output:
        config["out_dir"] + "/all_samples/pan_genome/samples.tsv"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        echo "sample\tgenome_fasta" > {output}
        echo "{input}" | tr ' ' '\n' | awk '{{split($0,a,"/"); print a[length(a)-2]"\t"$0}}' >> {output}
        """

rule iterative_map_to_pan_LQ:
    """
    Iteratively create LQ pan genome
    by mapping LQ samples to ref+HQ and
    adding novel sequences.
    """
    input:
       samples=config["out_dir"] + "/all_samples/pan_genome/samples.tsv",
       ref_genome=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genome.fasta",
       ref_gff=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes_no_alt_splicing.gff3",
       ref_proteins=config["out_dir"] + "/HQ_samples/HQ_pan/pan_proteins_no_alt_splicing.fasta"
    output:
       config["out_dir"] + "/all_samples/pan_genome/all_novel.fasta",
       config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta"
    params:
        map_to_pan_script=pipeline_dir + '/iterative_map_to_pan.py',
        min_len=config['min_length'],
        out_dir=config["out_dir"] + "/all_samples/pan_genome/",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/iterative_map_to_pan.yml'
    shell:
        """
         python {params.map_to_pan_script} {input.ref_genome} {input.ref_gff} {input.ref_proteins} {input.samples} {params.out_dir} --cpus {params.ppn} --min_len {params.min_len}
        """

rule prep_annotation_chunks:
    """
    Divide non-ref contigs into chunks for efficient parallel analysis
    """
    input:
        config["out_dir"] + "/all_samples/pan_genome/all_novel.fasta"
    output:
        config["out_dir"] + "/all_samples/non_ref/chunks/chunks.done"
    params:
        chunks=config['max_jobs'] - 1,
        out_dir=config["out_dir"] + "/all_samples/non_ref/chunks",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/faSplit.yml'
    shell:
        """
        faSplit sequence {input} {params.chunks} {params.out_dir}/chunk
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
        realpath {params.chunks_dir}/*.fa | awk '{{n=split($0,a,"/"); sub(".yml","",a[n]); print a[n]"\t"$0}}' >> {output}
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
        python {params.filter_script} {input.fasta} {input.fasta_map} {params.max_aed} {output}
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

rule simplify_nonref_protein_names:
    """
    Keep only PanGeneX part of names
    """
    input:
        config["out_dir"] + "/all_samples/annotation/non_redun_maker.proteins_filter_nodupl.fasta"
    output:
        config["out_dir"] + "/all_samples/annotation/non_redun_maker.proteins_filter_nodupl_simp.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        sed 's/>\(.*\)-R.*/>\\1/' {input} > {output}
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

rule create_pan_proteome:
    """
    Create pan genome proteins fasta
    """
    input:
        non_ref_proteins=config["out_dir"] + "/all_samples/annotation/non_redun_maker.proteins_filter_nodupl_simp.fasta",
        ref_plus_hq_proteins=config["out_dir"] + "/HQ_samples/HQ_pan/pan_proteins_no_alt_splicing.fasta"
    output:
        pan_proteome=config["out_dir"] + "/all_samples/pan_genome/pan_proteome.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        cat {input.ref_plus_hq_proteins} {input.non_ref_proteins} > {output.pan_proteome}
        """
        
rule create_pan_annotation:
    """
    Create pan genome annotation gff
    """
    input:
        non_ref_gff=config["out_dir"] + "/all_samples/annotation/non_redun_maker.genes_filter_nodupl.gff",
        ref_plus_hq_gff=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes_no_alt_splicing.gff3"
    output:
        pan_genes=config["out_dir"] + "/all_samples/pan_genome/pan_genes.gff"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        cat {input.ref_plus_hq_gff} {input.non_ref_gff} > {output.pan_genes}
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

rule map_HQ_genome_to_pan:
    """
    Map HQ genome sequence to pan
    genome in order to detect gene PAV
    """
    input:
        pan_genome=config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta",
        hq_genome=get_hq_sample_genome
    output:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sam"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/minimap2.yml'
    shell:
        """
        minimap2 -ax asm5 -t {params.ppn} {input.pan_genome} {input.hq_genome} -L > {output}
        """

rule HQ_sam_to_sorted_bam:
    """
    Convert SAM output for HQ samples to sorted BAM
    """
    input:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sam" 
    output:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sort.bam"
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

rule HQ_index_bam:
    input:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sort.bam"
    output:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sort.bam.bai"
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

rule pan_genes_gff_to_bed:
    """
    Convert pan mRNA features from GFF3
    to bed for easy handling
    """
    input:
        config["out_dir"] + "/all_samples/pan_genome/pan_genes.gff"
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_mRNA.bed"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        awk '$3 == "mRNA" {{split($9,a,";"); split(a[1],b,"="); print $1"\t"$4"\t"$5"\t"b[2]}}' {input} | sort -k1,1 -k2,2n > {output}
        """

rule calculate_HQ_mRNA_cov:
    """
    Calculate coverage of HQ genome
    sequences on pan genes.
    """
    input:
        bed=config["out_dir"] + "/all_samples/pan_genome/pan_mRNA.bed",
        bam=config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sort.bam"
    output:
        order=config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}.order",
        sorted_bed=config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/pan_mRNA.sort.bed",
        filtered_bam=config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.filter.sort.bam",
        cov=config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.bedCovHist"
    params:
        sort_script=utils_dir + "/custom_sort_bed.py",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/bedtools.yml'
    shell:
        """
        samtools view {input.bam} | cut -f3 | uniq > {output.order}
        python {params.sort_script} {input.bed} {output.order} > {output.sorted_bed}
        bedtools intersect -a {input.bam} -b {input.bed} -wa > {output.filtered_bam}
        bedtools coverage -a {output.sorted_bed} -b {output.filtered_bam} -hist -sorted > {output.cov}
        """

rule calculate_LQ_mRNA_cov:
    """
    Calculate coverage of LQ reads
    on pan genes.
    """
    input:
        bed=config["out_dir"] + "/all_samples/pan_genome/pan_mRNA.bed",
        bam=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sort.bam"
    output:
        order=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}.order",
        sorted_bed=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}//pan_mRNA.sort.bed",
        filtered_bam=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.filter.sort.bam",
        cov=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.bedCovHist"
    params:
        sort_script=utils_dir + "/custom_sort_bed.py",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/bedtools.yml'
    shell:
        """
        samtools view {input.bam} | cut -f3 | uniq > {output.order}
        python {params.sort_script} {input.bed} {output.order} > {output.sorted_bed}
        bedtools intersect -a {input.bam} -b {input.bed} -wa > {output.filtered_bam}
        bedtools coverage -a {output.sorted_bed} -b {output.filtered_bam} -hist -sorted > {output.cov}
        """

rule detect_HQ_PAV:
    """
    Determine per gene per sample PA
    based on sequence coverage
    """
    input:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.bedCovHist"
    output:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_PAV.tsv"
    params:
        detect_pav_script=pipeline_dir + '/detect_gene_PAV_from_bed_cov.py',
        min_frac_covered=config['HQ_min_cov'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.detect_pav_script} {input} 1 {params.min_frac_covered} {wildcards.sample} > {output}
        """

rule detect_LQ_PAV:
    """
    Determine per gene per sample PA
    based on sequence coverage
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.bedCovHist"
    output:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_PAV.tsv"
    params:
        detect_pav_script=pipeline_dir + '/detect_gene_PAV_from_bed_cov.py',
        min_depth=config['min_read_depth'],
        min_frac_covered=config['LQ_min_cov'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.detect_pav_script} {input} {params.min_depth} {params.min_frac_covered} {wildcards.sample}_{wildcards.ena_ref} > {output}
        """

rule create_ref_genes_list:
    """
    Create a file containing the
    reference mRNAs list
    """
    input:
        config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans_simp.gff.gene_to_mRNA'
    output:
        config["out_dir"] + "/all_samples/ref/ref_genes.list"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cut -f2 {input} > {output}
        """  

rule create_pan_PAV:
    """
    Combine HQ and LQ PAV tables
    and add ref column to create
    the final PAV matrix
    """
    input:
        lq=expand(config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_PAV.tsv", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        hq=expand(config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_PAV.tsv", sample=config['hq_info'].keys()),
        ref_genes=config["out_dir"] + "/all_samples/ref/ref_genes.list",
        names_sub=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans_simp.gff.gene_to_mRNA'
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv"
    params:
        create_PAV_matrix_script=pipeline_dir + '/create_PAV_matrix.py',
        ref_name=config['reference_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/pandas.yml'
    shell:
        """
        python {params.create_PAV_matrix_script} {input.hq} {input.lq} {params.ref_name} {input.ref_genes} {input.names_sub} {output}
        """
