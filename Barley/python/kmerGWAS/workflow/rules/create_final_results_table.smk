rule create_final_results_table:
    """
    Create a table with all
    significant k-mers, p-values
    direct mapping coordinates and
    projected coordinates on all LQ
    genomes
    """
    input:
        p_val_tsv = os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per_minus_logP.tsv'),
        direct_map_sam = os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_vs_%s.filter.sam' % ref_genome),
        projected_tsvs = expand(os.path.join(out_dir, 'all_samples', '{{phenotype}}', 'map_signif_kmers', 'pass_threshold_5per_vs_{genome}.filter.project_to_%s.tsv' % ref_genome), genome=LQ_genomes)
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per_mapping.tsv')
    params:
        create_table_script = os.path.join(scripts_dir, 'create_final_results_table.py')
    log:
        os.path.join(logs_dir, 'create_final_results_table', '{phenotype}.create_final_results_table.log')
    conda:
        os.path.join(envs_dir, 'project_mapping.yaml')
    resources:
        mem_gb = 8
    shell:
        """
        python {params.create_table_script} {input.p_val_tsv} {input.direct_map_sam} {input.projected_tsvs} {output}
        """
