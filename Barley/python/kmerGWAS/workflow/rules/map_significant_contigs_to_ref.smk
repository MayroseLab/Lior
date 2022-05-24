rule map_significant_contigs_to_ref:
    """
    Map significant contigs from LQ
    genomes to the reference genome 
    """
    input:
        contigs_fasta = os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_{genome}.signif_contigs.fasta'),
        ref_genome_fasta = os.path.join(out_dir, 'genomes', '%s.fasta' % ref_genome)
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_{genome}.signif_contigs_vs_%s.sam' % ref_genome)
    log:
        os.path.join(logs_dir, 'map_significant_contigs_to_ref', '{phenotype}_{genome}.map_significant_contigs_to_ref.log')
    conda:
        os.path.join(envs_dir, 'bwa.yaml')
    resources:
        mem_gb = 8
    shell:
        """
        bwa mem {input.ref_genome_fasta} {input.contigs_fasta} > {output}
        """
