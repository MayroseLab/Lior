rule create_introns_bed_flank_size_limit:
    """
    Create a BED file containing
    only intron with size <= PIR flank size
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.bed')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.flank_size_limit.bed')
    log:
        os.path.join(logs_dir, 'create_introns_bed_flank_size_limit', '{species}.create_introns_bed_flank_size_limit.log')
    params:
        flank_size = config['flank_size']
    shell:
        """
        awk '$3-$2 <= {params.flank_size}' {input} > {output}
        """
