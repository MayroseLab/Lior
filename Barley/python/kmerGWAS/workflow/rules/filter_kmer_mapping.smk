rule filter_sam:
    """
    Filter SAM of significant k-mers
    mapping to include only perfect
    mappings
    """
    input:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_vs_{genome}.sam')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_vs_{genome}.filter.sam')
    log:
        os.path.join(logs_dir, 'filter_kmer_mapping', '{phenotype}.{genome}.filter_sam.log')
    params:
        k = config['k']
    conda:
        os.path.join(envs_dir, 'samtools.yaml')
    resources:
        mem_gb = 16
    shell:
        """
        samtools view -h -e 'flag.unmap != 4 && qlen == {params.k} && [NM] && [NM] == 0' {input} > {output}
        """
