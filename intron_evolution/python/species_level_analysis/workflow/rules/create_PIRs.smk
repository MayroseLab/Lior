rule create_PIRs:
    """
    Create A BED with PIRs, randomly
    sampled from intergenic based
    on intron lengths
    """
    input:
        intergenic = os.path.join(out_dir, 'per_species', '{species}', 'intergenic.bed'),
        introns = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.bed')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'PIR.bed')
    log:
        os.path.join(logs_dir, 'create_PIRs', '{species}.create_PIRs.log')
    conda:
        os.path.join(envs_dir, 'biopython.yml')
    params:
        create_PIRs_script = os.path.join(scripts_dir, 'create_PIR_bed.py'),
        flank_size = config['flank_size']
    resources:
        mem_mb = 4000
    shell:
        """
        python {params.create_PIRs_script} {input.intergenic} {input.introns} {params.flank_size} {output}
        """
