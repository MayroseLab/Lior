rule project_kmers_to_ref:
    """
    Project k-mer mappings on contigs
    to reference coordinates
    """
    input:
        kmers = os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_vs_{genome}.filter.sam'),
        contigs = os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_{genome}.signif_contigs_vs_%s.filter.sam' % ref_genome)
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_vs_{genome}.filter.project_to_%s.tsv' % ref_genome)
    params:
        project_script = os.path.join(scripts_dir, 'project_kmer_mapping.py')
    log:
        os.path.join(logs_dir, 'project_kmers_to_ref', '{phenotype}_{genome}.project_kmers_to_ref.log')
    conda:
        os.path.join(envs_dir, 'project_mapping.yaml')
    resources:
        mem_gb = 8
    shell:
        """
        python {params.project_script} {input.kmers} {input.contigs} {wildcards.genome} {output}
        """
