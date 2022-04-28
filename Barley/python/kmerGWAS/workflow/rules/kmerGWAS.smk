rule kmerGWAS:
    """
    Perform GWAS per phenotype,
    using permutation-based threshold
    """
    input:
        kmer_matrix = os.path.join(out_dir, 'all_samples', 'kmers_PA_matrix.table'),
        kinship_matrix = os.path.join(out_dir, 'all_samples', 'kmers_PA_matrix.kinship'),
        phenotype_table = os.path.join(out_dir, 'all_samples', '{phenotype}', '{phenotype}.pheno')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'GWAS', 'kmers/output/phenotype_value.assoc.txt'),
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'GWAS', 'kmers/pass_threshold_5per')
    params:
        k = config['k'],
        kmerGWAS_script = os.path.join(config['voichek_code_dir'], 'kmers_gwas.py'),
        kmer_matrix_pref = os.path.join(out_dir, 'all_samples', 'kmers_PA_matrix'),
        out_dir = os.path.join(out_dir, 'all_samples', '{phenotype}', 'GWAS'),
        gemma_path = os.path.join(config['voichek_code_dir'], 'external_programs/gemma_0_96'),
        kmers_to_analyze = config['kmers_to_analyze']
    log:
        os.path.join(logs_dir, 'kmerGWAS', '{phenotype}.kmerGWAS.log')
    resources:
        mem_gb = 4
    threads:
        config['cpu']
    conda:
        os.path.join(envs_dir, 'kmerGWAS.yaml')
    shell:
        """
        rm -rf {params.out_dir}
        python {params.kmerGWAS_script} --pheno {input.phenotype_table} --kmers_table {params.kmer_matrix_pref} -l {params.k} -p {threads} --outdir {params.out_dir} -k {params.kmers_to_analyze}
        gzip -d {params.out_dir}/kmers/output/phenotype_value.assoc.txt.gz
        """
