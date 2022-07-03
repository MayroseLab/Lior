rule extract_singificant_kmers_fasta:
    """
    Extract sequences of all
    significant k-mers to fasta
    """
    input:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'GWAS', 'kmers/pass_threshold_5per')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per.fasta')
    log:
        os.path.join(logs_dir, 'extract_significant_kmers', '{phenotype}.extract_singificant_kmers_fasta.log')
    shell:
        """
        tail -n +2 {input} | awk '{{split($2,a,"_"); print ">kmer_"a[2]"\\n"a[1]}}' > {output}
        """

rule extract_singificant_kmers_list:
    """
    Extract sequences of all
    significant k-mers to simple list
    """
    input:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per.fasta')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per.list')
    log:
        os.path.join(logs_dir, 'extract_significant_kmers', '{phenotype}.extract_singificant_kmers_list.log')
    shell:
        """
        grep -v '>' {input} > {output} || true 
        """

rule extract_singificant_kmers_PAV:
    """
    Extract k-mer presence/absence
    for significant k-mes
    """
    input:
        lst = os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per.list'),
        matrix = os.path.join(out_dir, 'all_samples', 'kmers_PA_matrix.table')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per_PAV.tsv')
    params:
        filter_matrix_script = os.path.join(config['voichek_code_dir'], 'bin/filter_kmers'),
        matrix_pref = os.path.join(out_dir, 'all_samples', 'kmers_PA_matrix')
    log:
        os.path.join(logs_dir, 'extract_significant_kmers', '{phenotype}.extract_singificant_kmers_PAV.log')
    conda:
        os.path.join(envs_dir, 'kmerGWAS.yaml')
    resources:
        mem_gb = 4
    shell:
        """
        kmers=$(wc -l < {input.lst})
        if [ $kmers -gt 0 ]
        then
            {params.filter_matrix_script} -t {params.matrix_pref} -k {input.lst} -o {output}
        else
            names=$(tr '\\n' '\\t' < {params.matrix_pref}.names)
            echo -e "kmer\\t$names" > {output}
        fi
        """

rule  significant_kmers_minus_log_p:
    """
    Create a simplified table
    with all significant k-mers
    and their -log p-value
    """
    input:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'GWAS', 'kmers/pass_threshold_5per')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per_minus_logP.tsv')
    log:
        os.path.join(logs_dir, 'extract_significant_kmers', '{phenotype}.significant_kmers_minus_log_p.log')
    run:
        df = pd.read_csv(input[0], sep='\t')
        if df.empty:
            df = pd.DataFrame(columns=['k-mer_ID', 'k-mer_sequence', '-logP'])
        else:
            df[['k-mer_sequence', 'k-mer_ID']] = df['rs'].str.split('_', expand=True)
            df['k-mer_ID'] = 'kmer_' + df['k-mer_ID']
            df['-logP'] = -np.log10(df['p_lrt'])
            df = df[['k-mer_ID', 'k-mer_sequence', '-logP']].sort_values(by='-logP', ascending=False)
        df.to_csv(output[0], sep='\t', index=False)

rule significant_kmers_by_accession:
    """
    Create a TSV file with columns:
    Accession	# signif. k-mers	phenotype value
    sorted by # signif. k-mers (largest to smallest)
    """
    input:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per_PAV.tsv'),
        pheno_tsv = os.path.join(out_dir, 'all_samples', '{phenotype}', '{phenotype}.pheno')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per_by_acc.tsv')
    log:
        os.path.join(logs_dir, 'extract_significant_kmers', '{phenotype}.significant_kmers_by_accession.log')
    run:
        pav_df = pd.read_csv(input[0], sep='\t', index_col=0)
        pheno_df = pd.read_csv(input[1], sep='\t', index_col=0)
        kmers_per_acc = pd.DataFrame(pav_df.sum(axis=0))
        kmers_per_acc.columns = ['kmers_count']
        kmers_per_acc_with_pheno = pd.concat([kmers_per_acc, pheno_df], axis=1).dropna()
        kmers_per_acc_with_pheno.sort_values(by='kmers_count', ascending=False, inplace=True)
        kmers_per_acc_with_pheno.to_csv(output[0], sep='\t', index_label='Accession')
