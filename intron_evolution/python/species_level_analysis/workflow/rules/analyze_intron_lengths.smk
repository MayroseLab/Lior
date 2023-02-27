rule analyze_intron_lengths:
    """
    Analyze and extract stats
    about intron length from gff
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.gff3'),
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'intron_lengths.stats')
    log:
        os.path.join(logs_dir, 'analyze_intron_lengths', '{species}.analyze_intron_lengths.log')
    resources:
        mem_mb = 4000
    run:
        headers = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        gff_df = pd.read_csv(input[0], sep='\t', names=headers, comment='#')
        introns_df = gff_df.query('type == "intron"')
        introns_df['length'] = introns_df['end'] - introns_df['start'] + 1
        stats = introns_df['length'].describe()
        stats = stats[['min','max', 'mean', '50%', 'count']]
        stats.index = ['min', 'max', 'mean', 'median', 'count']
        mean_introns_per_trans = introns_df['attributes'].value_counts().mean()
        stats['mean_per_transcript'] = mean_introns_per_trans
        genes_df = gff_df.query('type == "mRNA"')
        genes_df['length'] = genes_df['end'] - genes_df['start'] + 1
        genes_total = genes_df['length'].sum()
        introns_total = introns_df['length'].sum()
        gene_intron_content = introns_total/genes_total*100
        stats['genes_intron_content'] = gene_intron_content
        stats['species'] = wildcards.species
        stats = stats[['species', 'min', 'max', 'mean', 'median', 'count', 'mean_per_transcript', 'genes_intron_content']]
        with open(output[0], 'w') as fo:
            print('\t'.join([str(s) for s in stats]), file=fo)
