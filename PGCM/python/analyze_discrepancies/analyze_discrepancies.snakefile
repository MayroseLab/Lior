"""
============================================================
This workflow performs all analyses
required for assigning causes to
discrepancies between a DN and a MTP
pan-genomes.
The main inputs are:
* the Panoramic output dirs for each
  pan genome
* a discrepancies TSV with columns:
  - gene
  - sample
  - gene type (ref/matched non-ref/unmatched non-ref)
  - type of discrepancy (DN+/MTP- or DN-/MTP+)
The output is the same table with an
additional column - discrepancy cause
============================================================
"""

import os
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))

LOGS_DIR = config['out_dir'] + "/logs"
CONDA_ENV_DIR = pipeline_dir + "/conda_env"

onstart:
    if not os.path.isdir(LOGS_DIR):
        os.makedirs(LOGS_DIR)

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

samples, ena_refs = [], []
for sample in os.listdir(config['dn_dir'] + '/per_sample/'):
    samples.append(sample)
    ena_ref = '_'.join(list(filter(lambda x: x.startswith('RG_assembly'), os.listdir(config['dn_dir'] + '/per_sample/' + sample)))[0].split('_')[2:])
    ena_refs.append(ena_ref)

wildcard_constraints:
    sample='[^_]+'

rule all:
    input:
        #config['out_dir'] + '/discrep_causes.tsv'
        expand(config['out_dir'] + '/MTP_trans_vs_assembly/{sample}_{ena_ref}.psl', zip, sample=samples, ena_ref=ena_refs),
        expand(config['out_dir'] + '/MTP_prot_vs_sample_proteins/{sample}_{ena_ref}.psl', zip, sample=samples, ena_ref=ena_refs),
        config['out_dir'] + '/DN_trans_vs_novel/' + 'DN_trans_vs_novel.psl',
        expand(config['out_dir'] + '/variant_calling/{sample}_{ena_ref}.filter.ann.bcf', zip, sample=samples, ena_ref=ena_refs),

rule search_mtp_transcripts_in_genomes:
    """
    search for all transcripts from the MTP
    pan-genome in each of the genome assemblies
    of the DN pan-genome
    """
    input:
        transcripts = config['mtp_dir'] + '/all_samples/pan_genome/pan_transcripts.fasta',
        assembly = config['dn_dir'] + '/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffolds.fasta'
    output:
        config['out_dir'] + '/MTP_trans_vs_assembly/{sample}_{ena_ref}.psl'
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/blat.yml'
    shell:
        """
        blat {input.assembly} {input.transcripts} {output}
        """

rule search_mtp_proteins_in_proteomes:
    """
    search for all proteins from the MTP
    pan-genome in each of the proteomes
    of the DN pan-genome
    """
    input:
        mtp_proteins = config['mtp_dir'] + '/all_samples/pan_genome/pan_proteome.fasta',
        sample_proteins = config['dn_dir'] + '/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter_nodupl.fasta'
    output:
        config['out_dir'] + '/MTP_prot_vs_sample_proteins/{sample}_{ena_ref}.psl'
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/blat.yml'
    shell:
        """
        blat {input.sample_proteins} {input.mtp_proteins} {output} -prot
        """

rule search_dn_transcripts_in_mtp_pan:
    """
    search for all transcripts from the DN
    pan-genome in the novel sequences of the
    MTP pan-genome
    """
    input:
        transcripts = config['dn_dir'] + '/all_samples/pan_genome/pan_transcripts.fasta',
        novel = config['mtp_dir'] + '/all_samples/pan_genome/all_novel.fasta'
    output:
        config['out_dir'] + '/DN_trans_vs_novel/' + 'DN_trans_vs_novel.psl'
    params:
        queue=config['queue'],
        priority=config['priority'],
        ppn=config['ppn'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/blat.yml'
    shell:
        """
        blat {input.novel} {input.transcripts} {output}
        """

rule index_reads_mapping:
    """
    Index the bam file that will
    be used for variant calling
    """
    input:
        config['mtp_dir'] + '/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.filter.sort.bam'
    output:
        config['mtp_dir'] + '/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.filter.sort.bam.bai'
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

rule index_mtp_pan:
    """
    Index the MTP pan-genome sequence fasta
    """
    input:
        config['mtp_dir'] + '/all_samples/pan_genome/pan_genome.fasta'
    output:
        config['mtp_dir'] + '/all_samples/pan_genome/pan_genome.fasta.fai'
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

rule call_variants:
    """
    Use reads mapping of each sample
    to call sequence variants
    """
    input:
        bam = config['mtp_dir'] + '/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.filter.sort.bam',
        bai = config['mtp_dir'] + '/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.filter.sort.bam.bai',
        pan_seq = config['mtp_dir'] + '/all_samples/pan_genome/pan_genome.fasta',
        pan_seq_fai = config['mtp_dir'] + '/all_samples/pan_genome/pan_genome.fasta.fai'
    output:
        config['out_dir'] + '/variant_calling/{sample}_{ena_ref}.bcf'
    params:
        queue=config['queue'],
        priority=config['priority'],
        ppn=config['ppn'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/bcftools.yml'
    shell:
        """
        bcftools mpileup -Ou -f {input.pan_seq} {input.bam} | bcftools call -mv -Ob -o {output} --threads {params.ppn}
        """

rule filter_variants:
    """
    Discard low quality variant calls
    """
    input:
        config['out_dir'] + '/variant_calling/{sample}_{ena_ref}.bcf'
    output:
        config['out_dir'] + '/variant_calling/{sample}_{ena_ref}.filter.bcf'
    params:
        min_quality=config['min_qualuty'],
        min_depth=config['min_depth'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/bcftools.yml'
    shell:
        """
        bcftools filter --exclude '%QUAL<{params.min_quality} || DP<{params.min_depth}' {input} -o {output}
        """

rule build_snpeff_db:
    """
    Build a SNPEff DB of the MTP pan-genome
    """
    input:
        gff = config['mtp_dir'] + '/all_samples/pan_genome/pan_genes.gff',
        genome = config['mtp_dir'] + '/all_samples/pan_genome/pan_genome.fasta',
        prot = config['mtp_dir'] + '/all_samples/pan_genome/pan_proteome.fasta'
    output:
        config['out_dir'] + '/SNPEff_db.done'
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snpeff.yml'
    shell:
        """
        # find snpeff.config
        env=`grep -l snpeff ./.snakemake/conda/*.yaml | xargs ls -1tr | tail -1 | xargs basename | sed 's/\.yaml//'`
        config="./.snakemake/conda/$env/share/snpeff-5.0-0/snpEff.config"
        # make sure file exists
        if [ ! -f "$config" ]; then
            echo "Can't find SNPEff config $config"
            exit 1
        fi
        # If pan-genome DB is not already there, add to config
        if ! grep -qw pan-genome.genome $config; then
            echo 'pan-genome.genome : pan-genome' >> $config
        fi
        # get pan-genome gff, proteins and sequences data
        dataDir="./.snakemake/conda/$env/share/snpeff-5.0-0/data/pan-genome"
        mkdir -p $dataDir
        ln -f {input.gff} $dataDir/genes.gff
        ln -f {input.genome} $dataDir/sequences.fa
        ln -f  {input.prot} $dataDir/protein.fa
        # Build DB
        snpEff build -gff3 -v pan-genome
        touch {output}
        """

rule annotate_variants:
    """
    Use SNPEff to annotate variants
    """
    input:
        bcf = config['out_dir'] + '/variant_calling/{sample}_{ena_ref}.filter.bcf',
        db = config['out_dir'] + '/SNPEff_db.done'
    output:
        config['out_dir'] + '/variant_calling/{sample}_{ena_ref}.filter.ann.bcf'
    params:
        queue=config['queue'],
        priority=config['priority'],
        ppn=2,
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snpeff.yml'
    shell:
        """
        snpEff pan-genome {input.bcf} > {output}
        """
 
