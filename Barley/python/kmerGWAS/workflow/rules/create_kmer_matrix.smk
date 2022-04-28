rule create_kmers_path_list:
    """
    Create list of combined k-mer count
    files per accession
    """
    input:
        expand(os.path.join(out_dir, 'per_sample', '{sample}','{sample}_kmers_with_strand'), sample=samples)
    output:
        os.path.join(out_dir, 'all_samples', 'kmers_list_paths.txt')
    log:
        os.path.join(logs_dir, 'create_kmer_matrix', 'create_kmers_path_list.log')
    run:
        with open(output[0],'w') as fo:
            for p in input:
                sample = os.path.basename(p).replace('_kmers_with_strand','')
                print("%s\t%s" %(p, sample), file=fo)

rule create_kmers_list:
    """
    Create a filtered list of k-mers to be
    included in the k-mers matrix.
    """
    input:
        path_list = os.path.join(out_dir, 'all_samples', 'kmers_list_paths.txt'),
        files = expand(os.path.join(out_dir, 'per_sample', '{sample}','{sample}_kmers_with_strand'), sample=samples)
    output:
        temp(os.path.join(out_dir, 'all_samples', 'kmers_to_use'))
    params:
        k = config['k'],
        min_accessions = config['min_accessions'],
        min_perc_acc = config['min_perc_acc'],
        list_kmers_script = os.path.join(config['voichek_code_dir'], 'bin/list_kmers_found_in_multiple_samples')
    log:
        os.path.join(logs_dir, 'create_kmer_matrix', 'create_kmers_list.log')
    resources:
        mem_gb = config['mem_gb']
    conda:
        os.path.join(envs_dir, 'gcc.yaml')
    shell:
        """
        {params.list_kmers_script} -l {input.path_list} -k {params.k} --mac {params.min_accessions} -p {params.min_perc_acc} -o {output}
        """

rule create_kmers_matrix:
    """
    Create filtered matrix of k-mer
    presence-absence across accessions
    """
    input:
        paths_list = os.path.join(out_dir, 'all_samples', 'kmers_list_paths.txt'),
        kmers_list = os.path.join(out_dir, 'all_samples', 'kmers_to_use')
    output:
        os.path.join(out_dir, 'all_samples', 'kmers_PA_matrix.table')
    params:
        k = config['k'],
        build_kmers_table_script = os.path.join(config['voichek_code_dir'], 'bin/build_kmers_table'),
        out_pref = os.path.join(out_dir, 'all_samples', 'kmers_PA_matrix')
    log:
        os.path.join(logs_dir, 'create_kmer_matrix', 'create_kmers_matrix.log')
    resources:
        mem_gb = config['mem_gb']
    conda:
        os.path.join(envs_dir, 'gcc.yaml')
    shell:
        """
        {params.build_kmers_table_script} -l {input.paths_list} -k {params.k} -a {input.kmers_list} -o {params.out_pref}
        """
