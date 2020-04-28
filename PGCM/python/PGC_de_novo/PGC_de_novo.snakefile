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

localrules: all, prep_liftover_chunks_tsv, prep_annotation_chunks_tsv, prep_liftover_yaml, prep_annotation_yaml, require_evidence

rule all:
    input:
        pav=config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv",
        cnv=config["out_dir"] + "/all_samples/pan_genome/pan_CNV.tsv",
        prot=config["out_dir"] + "/all_samples/pan_genome/pan_proteins.fasta"

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
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs.fasta"
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

rule index_contigs_fasta:
    """
    Index contigs fasta. For later use by RaGOO.
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta.fai"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools faidx {input}
        """

rule assembly_quast:
    """
    Run QUAST on assembly to get assembly stats and QA
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

rule ref_guided_assembly:
    """
    Assemble contigs into pseudomolecules
    by mapping to the reference genome and
    using reference-guided assembly
    """
    input:
        contigs=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta",
        ref_genome=config['ref_genome'],
    output:
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/ragoo.fasta",
        directory(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/orderings")
    params:
        ragoo_script=config['ragoo_script'],
        gap_size=config['gap_size'],
        out_dir=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}",
        queue=config['queue'],
        priority=config['priority'],
        ppn=config['ppn'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/RaGOO.yml'
    shell:
        """
        cd {params.out_dir}
        ln {input.contigs} contigs.fasta
        ln {input.ref_genome} ref.fasta
        python {params.ragoo_script} contigs.fasta ref.fasta -t {params.ppn} -g {params.gap_size}
        """

rule assembly_busco:
    """
    Run BUSCO on assembly
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/ragoo.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/run_BUSCO/short_summary_BUSCO.txt"
    params:
        assembly_dir=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output",
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

rule prep_chunks:
    """
    Divide assembly into chunks for efficient parallel analysis
    """
    input:
        fasta=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/ragoo.fasta",
    output:
        config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/chunks.lft"
    params:
        n_chunks=config['max_jobs']//len(config['samples_info']) - 1,
        out_pref=config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/chunk",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/faSplit.yml'
    shell:
        """
        chunkSize=`expr $(grep -v '>' {input} | tr -d '\n' | wc | awk '{{print $3}}') / {params.n_chunks}`
        faSplit gap {input} $chunkSize {params.out_pref} -noGapDrops -minGapSize=10 -lift={output}
        """ 

rule prep_liftover_chunks_tsv:
    """
    Prepare TSV config for liftover run
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/chunks.lft"
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
        realpath {params.chunks_dir}/*.fa | awk '{{n=split($0,a,"/"); print a[n]"\t"$0}}' >> {output}
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
        echo "sample: {wildcards.sample}" >> {output}
        echo "logs_dir: {params.logs_dir}" >> {output}
        echo config_kv_pairs: est={params.liftover_transcripts},{params.additional_transcripts} protein={params.proteins} rmlib={params.repeats_library} maker_gff={input.liftover_gff} augustus_species={params.augustus_species} >> {output}
        """

rule maker_annotation:
    """
    Run MAKER using liftover results + evidence
    to obtain a full genome annotation
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/annotation.yml"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.fasta",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.all.gff"
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

rule rename_genes:
    """
    Assign genes short, unique names (gff and fasta).
    Names consist of the genome name and a unique ID.
    """
    input:
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.fasta",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.gff"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.raw_names.fasta",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.raw_names.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/gff.map",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/fasta.map"
    params:
        out_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/",
        create_fasta_map_script=os.path.join(utils_dir,"create_fasta_map.py"),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        maker_map_ids --prefix {wildcards.sample}_ --justify 1 --iterate 1 {input.gff} > {params.out_dir}/gff.map
        cp {input.gff} {params.out_dir}/maker.genes.raw_names.gff
        map_gff_ids {params.out_dir}/gff.map {input.gff}
        python {params.create_fasta_map_script} {input.gff} > {params.out_dir}/fasta.map
        cp {input.fasta} {params.out_dir}/maker.proteins.raw_names.fasta
        map_fasta_ids {params.out_dir}/fasta.map {input.fasta}
        """

rule make_chunks_bed:
    """
    Create a bed file containing chunks borders.
    Useful as part of annotation evidence collection.
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/chunks.lft"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/chunks.bed"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        awk '{{print $4"\t"$1"\t"$1+$3"\t"$2}}' {input} > {output}
        """

rule convert_chunks_to_chromosomes:
    """
    Convert gff coordinates and sequence names
    to transform from chunks to chromosomes
    """
    input:
        genes=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.filter.gff",
        all_=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.all.gff",
        bed=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/chunks.bed"
    output:
        genes_out=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.filter.chr.gff",
        all_out=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.all.chr.gff"
    params:
        convert_script=os.path.join(utils_dir,"transform_gff_coordinates.py"),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.convert_script} {input.genes} {input.bed} {output.genes_out}
        python {params.convert_script} {input.all_} {input.bed} {output.all_out}
        """

rule make_evidence_gffs:
    """
    Extract evidence features from MAKER gff.
    Useful as part of annotation evidence collection.
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.all.chr.gff"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.augustus.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.blastn.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.blastx.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.est2genome.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.pred_gff:maker.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.protein2genome.gff"
    params:
        split_gff_script=os.path.join(utils_dir,"split_gff_by_source.py"),
        out_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.split_gff_script} {input} augustus,blastn,blastx,est2genome,pred_gff:maker,protein2genome {params.out_dir} 
        """

rule make_contigs_bed:
    """
    Create a bed file with original contig borders.
    Useful as part of annotation evidence collection.
    """
    input:
        ord_dir=directory(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/orderings"),
        faidx=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta.fai"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/contigs.bed"
    params:
        ragoo_contigs_script=os.path.join(pipeline_dir,"ragoo_ordering_to_bed.py"),
        gap_size=config['gap_size'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.ragoo_contigs_script} {input.ord_dir} {input.faidx} {params.gap_size} {output}
        """

rule index_evidence_gff:
    """
    Sort, compress and index (tabix)
    evidence gffs for easy display in IGV
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/contigs.bed",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.augustus.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.blastn.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.blastx.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.est2genome.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.pred_gff:maker.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.protein2genome.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/chunks.bed",
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/QUAST/report.html",
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/run_BUSCO/short_summary_BUSCO.txt"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/done"
    params:
        evidence_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/index_gff.yml'
    shell:
        """
        cd {params.evidence_dir}
        for x in `ls -1 *.gff`; do srt=`echo $x | sed 's/\.gff/\.sort\.gff/'`; bedtools sort -i $x > $srt; bgzip $srt; tabix $srt.gz; done
        touch {output}
        """

rule filter_annotation:
    """
    Remove unreliable genes from gff.
    These are genes that do NOT come from
    lift-over AND have AED > X (set by user)
    """
    input:
        gff_map=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/gff.map",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.gff"
    output:
        lst=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.filter.list",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.filter.gff"
    params:
        max_aed=config['max_aed'],
        filter_gff_script=utils_dir + '/filter_gff_by_id_list.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        awk '{{split($9,a,";"); split(a[1],b,"="); split(a[2],c,"="); split(a[5],d,"=")}} $3 == "mRNA" && (a[1] ~ /pred_gff/ || d[2] <= {params.max_aed}) {{print(b[2]"\\n"c[2])}}' {input.gff} > {output.lst}
        python {params.filter_gff_script} {input.gff} {output.lst} {output.gff}
        """

rule filter_proteins:
    """
    Filter proteins fasta according to filtered gff  
    """
    input:
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.filter.gff",
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.fasta",
        fasta_map=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/fasta.map"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter.fasta"
    params:
        filter_fasta_script=utils_dir + '/filter_fasta_by_gff.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.filter_fasta_script} {input.gff} {input.fasta} {output} mRNA ID
        """

rule prevent_duplicate_names:
    """
    If for any reason the filtered
    proteins contain duplicate names,
    prevent this by renaming. This
    helps prevent problems in next steps
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter_nodupl.fasta"
    params:
        dupl_script=utils_dir + '/prevent_duplicate_names.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.dupl_script} {input} {output}
        """

rule annotation_busco:
    """
    Run BUSCO on filtered annotation proteins
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter_nodupl.fasta"
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

rule prep_for_orthofinder:
    """
    Prepare orthofinder input - simplify
    fasta record names and put all fasta
    files into one dir with file names
    matching genome names.
    """
    input:
        ev=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/done",
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter_nodupl.fasta"
    output:
        config["out_dir"] + "/all_samples/orthofinder/{sample}_{ena_ref}.fasta"
    params:
        of_dir=config["out_dir"] + "/all_samples/orthofinder",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        sed 's/ protein .*//' {input.fasta} > {output}
        """

rule remove_ref_alt_splicing:
    """
    In case the reference annotation
    contains genes with multiple mRNAs,
    only keep the longest transcript.
    """
    input:
        config['ref_annotation']
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

rule get_ref_proteins:
    """
    Filter reference proteins according
    to filtered gff and put the new file
    in orthofinder dir.
    """
    input:
        fasta=config['ref_proteins'],
        gff=config["out_dir"] + "/all_samples/ref/" + config['ref_name'] + '_longest_trans.gff'
    output:
        config["out_dir"] + "/all_samples/orthofinder/" + config['ref_name'] + '.fasta'
    params:
        filter_fasta_script=utils_dir + '/filter_fasta_by_gff.py',
        name_attribute=config['name_attribute'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.filter_fasta_script} {input.gff} {input.fasta} {output} mRNA {params.name_attribute}
        """

rule orthofinder:
    """
    Run OrthoFinder2 on all proteins
    from all annotated genomes to get
    initial orthogroups
    """
    input:
        expand(config["out_dir"] + "/all_samples/orthofinder/{sample}_{ena_ref}.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        config["out_dir"] + "/all_samples/orthofinder/" + config['ref_name'] + '.fasta'
    output:
        config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups/Orthogroups.tsv"
    params:
        orthofinder_dir=config["out_dir"] + "/all_samples/orthofinder",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/orthofinder2.yml'
    shell:
        """
        if [ -d {params.orthofinder_dir}/OrthoFinder/Results_orthofinder/  ]
        then
            rm -rf {params.orthofinder_dir}/OrthoFinder/Results_orthofinder/
        fi
        orthofinder -t {params.ppn} -a {params.ppn} -S diamond -n orthofinder -f {params.orthofinder_dir}
        """

rule break_orthogroups_MWOP:
    """
    Use orthofinder's ouput and further
    break orthogroups into smaller
    clusters representing single genes.
    This is done using the Maximum Weight
    Orthogonal Partitions (MWOP) algorithm.
    """
    input:
        config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups/Orthogroups.tsv"
    output:
        config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups_break_MWOP.tsv"
    params:
        orthofinder_dir=config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder",
        mwop_script=os.path.join(pipeline_dir,"break_OrthoFinder_clusters.py"),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/break_orthogroups.yml'
    shell:
        """
        python {params.mwop_script} {params.orthofinder_dir} bitscore {output}
        """

rule create_PAV_matrix:
    """
    Create the final PAV and CNV matrices
    """
    input:
        config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups_break_MWOP.tsv"
    output:
        pav=config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv",
        cnv=config["out_dir"] + "/all_samples/pan_genome/pan_CNV.tsv",
        mapping=config["out_dir"] + "/all_samples/pan_genome/OG_to_gene_names.tsv"
    params:
        create_pav_mat_script=os.path.join(pipeline_dir,"create_PAV_matrix.py"),
        ref_name=config['ref_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/break_orthogroups.yml'
    shell:
        """
        python {params.create_pav_mat_script} {input} {params.ref_name} {output.pav} {output.cnv} {output.mapping}
        """

rule create_pan_proteins_fasta:
    """
    Create a fasta file with one
    representative protein sequence
    per pan gene.
    """
    input:
        mapping=config["out_dir"] + "/all_samples/pan_genome/OG_to_gene_names.tsv",
        mwop=config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups_break_MWOP.tsv"
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_proteins.fasta"
    params:
        create_pan_prot_fasta_script=os.path.join(pipeline_dir,"create_pan_proteins_fasta.py"),
        og_seq_dir=config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroup_Sequences",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        python {params.create_pan_prot_fasta_script} {params.og_seq_dir} {input.mapping} {input.mwop} {output}
        """
