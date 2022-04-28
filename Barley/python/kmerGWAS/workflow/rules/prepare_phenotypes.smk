out_dir = config['out_dir']

rule prepare_phenotype_table:
    """
    Prepare phenotype tables for GWASש
    """
    input:
        config['phenotypes_tsv']
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', '{phenotype}.pheno')
    log:
        os.path.join(logs_dir, '{phenotype}_prepare_phenotype_table', 'prepare_phenotypes.log')
    run:
        pheno = wildcards.phenotype
        one_pheno_df = pd.DataFrame(phenotypes_df[pheno])
        one_pheno_df = one_pheno_df.dropna()
        one_pheno_df.columns = ['phenotype_value']
        one_pheno_df.index.name = 'accession_id'
        pheno_dir = os.path.join(out_dir, 'all_samples', pheno)
        if not os.path.isdir(pheno_dir):
            os.mkdir(pheno_dir)
        out_file = os.path.join(pheno_dir, '%s.pheno' % pheno)
        one_pheno_df.to_csv(out_file, sep='\t', index=True)

rule plot_phenotype_histogram:
    """
    Create a HTML with histogram
    of phenotype values
    """
    input:
        os.path.join(out_dir, 'all_samples', '{phenotype}', '{phenotype}.pheno')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', '{phenotype}_histogram.html')
    log:
        os.path.join(logs_dir, 'prepare_phenotypes', '{phenotype}_plot_phenotype_histogram.log')
    params:
        histogram_script = os.path.join(scripts_dir, 'phenotype_histogram.py')
    conda:
        os.path.join(envs_dir,'plotly.yaml')
    shell:
        """
        python {params.histogram_script} {input} {wildcards.phenotype} {output}
        """
