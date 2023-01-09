rule repeat_content_intron_length_model:
    """
    Linear model: intron length ~ repeat content
    """
    input:
        rep_per_intron = os.path.join(out_dir, 'per_species', '{species}', 'repeats_per_intron.tsv'),
        introns_bed = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.bed')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'repeat_content_intron_length_model.tsv')
    log:
        os.path.join(logs_dir, 'repeat_content_intron_length_model', '{species}.repeat_content_intron_length_model.log')
    run:
        rep_per_intron_df = pd.read_csv(input['rep_per_intron'], sep='\t', index_col=0)
        introns_bed_df = pd.read_csv(input['introns_bed'], sep='\t', names=['seqid', 'start', 'end'], index_col=3)
        rep_per_intron_df = pd.concat([rep_per_intron_df, introns_bed_df], axis=1)
        rep_per_intron_df['length'] = rep_per_intron_df['end'] - rep_per_intron_df['start']
        if rep_per_intron_df.shape[0] < 30:
            slope, intercept, r_squared_value, p_value, std_err = ['nan']*5
        else:
            slope, intercept, r_value, p_value, std_err = stats.linregress(rep_per_intron_df['% covered by repeats'], rep_per_intron_df['length'])
            r_squared_value = r_value**2
        values = [wildcards.species] + [str(v) for v in [slope, intercept, r_squared_value, p_value, std_err]]
        with open(output[0], 'w') as fo:
            print('\t'.join(values), file=fo)
