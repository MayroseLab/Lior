rule create_final_report:
    """
    Create an HTML with various
    figures to summarize the
    analysis
    """
    input:
        pheno_tsv = os.path.join(out_dir, 'all_samples', '{phenotype}', '{phenotype}.pheno'),
        signif_kmers_by_acc_tsv = os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per_by_acc.tsv'),
        kmer_mapping_tsv = os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per_mapping.tsv'),
        chr_len_tsv = os.path.join(out_dir, 'genomes', '%s.chr.len' % ref_genome),
        ref_gff = config['reference_annotation_gff'],
        genomes_tsv = config['genomes_tsv']
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'report.html')
    params:
        create_report_script = os.path.join(scripts_dir, 'create_phenotype_report.py'),
        ref_name = ref_genome,
        contig_mapping_dir = os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers')
    log:
        os.path.join(logs_dir, 'create_final_report', '{phenotype}.create_final_report.log')
    conda:
        os.path.join(envs_dir, 'plotly.yaml')
    resources:
        mem_gb = 8
    shell:
        """
        python {params.create_report_script} {wildcards.phenotype} {input.pheno_tsv} {input.signif_kmers_by_acc_tsv} {input.kmer_mapping_tsv} {params.ref_name} {input.chr_len_tsv} {input.ref_gff} {params.contig_mapping_dir} {input.genomes_tsv} 1000000 5 {output}
        """
