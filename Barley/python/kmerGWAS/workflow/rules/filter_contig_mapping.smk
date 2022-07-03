rule filter_contig_mapping:
    """
    Remove unmapped and multiple-mappings
    """
    input:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_{genome}.signif_contigs_vs_%s.sam' % ref_genome)
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_{genome}.signif_contigs_vs_%s.filter.sam' % ref_genome)
    params:
        filter_script = os.path.join(scripts_dir, 'filter_contig_mapping.py')
    log:
        os.path.join(logs_dir, 'filter_contig_mapping', '{phenotype}_{genome}.filter_contig_mapping.log')
    conda:
        os.path.join(envs_dir, 'plotly.yaml')
    resources:
        mem_gb = 8
    shell:
        """
        python {params.filter_script} {input} 30 {output}
        """
