rule collect_repeat_stats:
    """
    Collect repeat stats from
    all species into one table
    """
    input:
        stats = expand(os.path.join(out_dir, 'per_species', '{species}', 'repeats.stats'), species=species_list),
        model = expand(os.path.join(out_dir, 'per_species', '{species}', 'repeat_content_intron_length_model.tsv'), species=species_list)
    output:
        stats = os.path.join(out_dir, 'all_species', 'repeats.stats'),
        model = os.path.join(out_dir, 'all_species', 'repeats_intron_length_model.stats')
    log:
        os.path.join(logs_dir, 'collect_repeat_stats', 'all_species.collect_repeat_stats.log')
    run:
        shell("cat {input.stats} > {output.stats}")
        headers = ['Repeat type', 'Assembly size', 'Total repeas', 'Total introns', 'Total intergenic', 'Percent of genome covered by repeats', 'Percent of introns covered by repeats', 'Percent of intergenic covered by repeats', 'Mean repeat length in introns', 'Mean repeat length in introns shorter than flank size', 'Mean repeat length in intergenic']
        rep_df = pd.read_csv(output['stats'], sep='\t', names=headers, index_col=0)
        rep_df['group'] = rep_df.index.map(species_df['group'])
        rep_df.to_csv(output['stats'], sep='\t', float_format=lambda x: round(x,2))

        shell("cat {input.model} > {output.model}")
        headers = ['Slope', 'Intercept', 'R^2', 'P value', 'STD Err']
        model_df = pd.read_csv(output['model'], sep='\t', names=headers, index_col=0)
        model_df = pd.concat([species_df, model_df], axis=1)
        model_df.to_csv(output['model'], sep='\t', float_format=lambda x: round(x,2))

