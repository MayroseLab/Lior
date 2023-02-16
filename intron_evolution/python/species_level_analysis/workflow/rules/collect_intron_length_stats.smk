rule collect_intron_length_stats:
    """
    Collect intron lengths stats from
    all species into one table
    """
    input:
        intron_length = expand(os.path.join(out_dir, 'per_species', '{species}', 'intron_lengths.stats'), species=species_list),
        genome_size = expand(os.path.join(out_dir, 'per_species', '{species}', 'genome.size'), species=species_list)
    output:
        os.path.join(out_dir, 'all_species', 'intron_lengths.stats')
    log:
        os.path.join(logs_dir, 'collect_intron_length_stats', 'all_species.collect_intron_length_stats.log')
    run:
        shell("cat {input.intron_length} > {output}.il")
        shell("cat {input.genome_size} > {output}.gs")
        headers = ['Min', 'Max', 'Mean', 'Median', 'Total number of introns', 'Mean number of introns per transcript']
        len_df = pd.read_csv(output[0]+'.il', sep='\t', names=headers, index_col=0)
        gs_df = pd.read_csv(output[0]+'.gs', sep='\t', names=['genome_size'], index_col=0)
        len_df = pd.concat([species_df, len_df, gs_df], axis=1)
        len_df.to_csv(output[0], sep='\t', float_format=lambda x: round(x,2))
