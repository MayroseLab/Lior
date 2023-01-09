rule introns_repeat_coverage:
    """
    Calculate the coverage of
    repeats on introns
    """
    input:
        introns = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.bed'),
        repeats = os.path.join(out_dir, 'per_species', '{species}', 'repeats', 'all.merge.bed')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'repeats_on_introns.cov.hist')
    log:
        os.path.join(logs_dir, 'introns_repeat_coverage', '{species}.introns_repeat_coverage.log')
    conda:
        os.path.join(envs_dir, 'bedtools.yml')
    resources:
        mem_mb = 4000
    shell:
        """
        bedtools coverage -a {input.introns} -b {input.repeats} -hist > {output}
        """
