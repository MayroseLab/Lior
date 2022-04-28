rule create_kinship_matrix:
    """
    Create kinship matrix from k-mers
    to control for population structure
    """
    input:
        os.path.join(out_dir, 'all_samples', 'kmers_PA_matrix.table')
    output:
        os.path.join(out_dir, 'all_samples', 'kmers_PA_matrix.kinship')
    params:
        k = config['k'],
        matrix_pref = os.path.join(out_dir, 'all_samples', 'kmers_PA_matrix'),
        min_maf = config['min_maf_for_kinship'],
        build_kinship_matrix_script = os.path.join(config['voichek_code_dir'], 'bin/emma_kinship_kmers')
    log:
        os.path.join(logs_dir, 'create_kinship_matrix', 'create_kinship_matrix.log')
    resources:
        mem_gb = config['mem_gb']
    conda:
        os.path.join(envs_dir, 'gcc.yaml')
    shell:
        """
        {params.build_kinship_matrix_script} -t {params.matrix_pref} -k {params.k} --maf {params.min_maf} > {output}
        """

