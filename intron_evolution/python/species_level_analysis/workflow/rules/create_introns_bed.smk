rule create_introns_bed:
    """
    Create BED file with
    intron features, including
    intron IDs
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.gff3')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.bed')
    log:
        os.path.join(logs_dir, 'create_introns_bed', '{species}.create_introns_bed.log')
    run:
        headers = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        gff_df = pd.read_csv(input[0], sep='\t', names=headers, comment='#')
        introns_df = gff_df.query('type == "intron"')
        introns_df['transcript'] = introns_df['attributes'].str.replace('Parent=','')
        introns_plus_df = introns_df.query('strand == "+"')
        introns_minus_df = introns_df.query('strand == "-"')
        introns_plus_df['intron_n'] = introns_plus_df.groupby('transcript').cumcount()
        introns_minus_df['intron_n'] = introns_minus_df.groupby('transcript').cumcount(ascending=False)
        introns_df = pd.concat([introns_plus_df, introns_minus_df])
        introns_df['intron_id'] = introns_df['transcript'] + '_' + introns_df['intron_n'].astype(str)
        introns_bed_df = introns_df[['seqid', 'start','end', 'intron_id']]
        introns_bed_df['start'] = introns_bed_df['start'] - 1
        introns_bed_df.to_csv(output[0], sep='\t', header=False, index=False)
