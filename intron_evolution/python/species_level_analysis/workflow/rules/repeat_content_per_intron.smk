rule repeat_content_per_intron:
    """
    Calculate the % covered
    by repeats for each intron.
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'repeats_on_introns.cov.hist')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'repeats_per_intron.tsv')
    log:
        os.path.join(logs_dir, 'repeat_content_per_intron', '{species}.repeat_content_per_intron.log')
    resources:
        mem_mb = 4000
    run:
        header = ['seqid', 'start', 'end', 'intron_id', 'repeat_cov', 'cov_length', 'intron_length', 'cov_frac']
        cov_hist_df = pd.read_csv(input[0], sep='\t', names=header)
        cov_hist_df = cov_hist_df.query('seqid != "all"')

        def frac_cov(df):
            cov = df.query('repeat_cov == 1')
            if cov.empty:
                return pd.Series([0], index=[1])
            return pd.Series([cov.iloc[0]['cov_frac']*100], index=[1])

        repeat_cov_per_intron = cov_hist_df.groupby('intron_id').apply(frac_cov)
        repeat_cov_per_intron.columns = ['% covered by repeats']
        repeat_cov_per_intron.index.name = 'Intron ID'
        repeat_cov_per_intron.to_csv(output[0], sep='\t')
