rule create_bedtools_genome:
    """
    Create a genome file to
    be used by bedtools
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'genome.sm.fasta')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'genome.sm.genome')
    log:
        os.path.join(logs_dir, 'create_bedtools_genome', '{species}.create_bedtools_genome.log')
    conda:
        os.path.join(envs_dir, 'samtools.yml')
    resources:
        mem_mb = 4000
    shell:
        """
        samtools faidx {input}
	cut -f1,2 {input}.fai > {output}
        """
