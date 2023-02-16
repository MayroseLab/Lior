rule calculate_genome_size:
    """
    Calculate and write out
    genome (assembly) size
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'genome.sm.genome')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'genome.size')
    log:
        os.path.join(logs_dir, 'calculate_genome_size', '{species}.calculate_genome_size.log')
    shell:
        """
        gs=$(awk '{{SUM+=$2}}END{{print SUM}}' {input})
        echo -e "{wildcards.species}\t$gs" > {output}
        """
