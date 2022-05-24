rule select_significant_contigs:
    """
    Create a list of contigs from LQ
    genomes containing at least X
    significant k-mers
    """
    input:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_vs_{genome}.filter.sam')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_{genome}.signif_contigs.list')
    params:
        min_kmers = config['min_kmers']
    log:
        os.path.join(logs_dir, 'extract_significant_contigs', '{phenotype}_{genome}.select_significant_contigs.log')
    conda:
        os.path.join(envs_dir, 'samtools.yaml')
    resources:
        mem_gb = 16
    shell:
        """
        samtools view {input} | cut -f3 | sort | uniq -c | awk '$1 >= {params.min_kmers} {{print $2}}' > {output}
        """

rule extract_significant_contigs:
    """
    Extract FASTA sequences of contigs
    containing at least X significant k-mers
    """
    input:
        kmer_list = os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_{genome}.signif_contigs.list'),
        genome_fasta = get_genome
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_{genome}.signif_contigs.fasta')
    log:
        os.path.join(logs_dir, 'extract_significant_contigs', '{phenotype}_{genome}.extract_significant_contigs.log')
    conda:
        os.path.join(envs_dir, 'seqtk.yaml')
    resources:
        mem_gb = 16
    shell:
        """
        seqtk subseq {input.genome_fasta} {input.kmer_list} > {output}
        """
