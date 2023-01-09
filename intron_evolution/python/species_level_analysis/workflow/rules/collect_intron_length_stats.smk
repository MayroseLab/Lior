rule collect_intron_length_stats:
    """
    Collect intron lengths stats from
    all species into one table
    """
    input:
        expand(os.path.join(out_dir, 'per_species', '{species}', 'intron_lengths.stats'), species=species_list)
    output:
        os.path.join(out_dir, 'all_species', 'intron_lengths.stats')
    log:
        os.path.join(logs_dir, 'collect_intron_length_stats', 'all_species.collect_intron_length_stats.log')
    run:
        shell("cat {input} > {output}")
        headers = ['Min', 'Max', 'Mean', 'Median', 'Total number of introns', 'Mean number of introns per transcript']
        len_df = pd.read_csv(output[0], sep='\t', names=headers, index_col=0)
        len_df = pd.concat([species_df, len_df], axis=1)
        len_df.to_csv(output[0], sep='\t', float_format=lambda x: round(x,2))
