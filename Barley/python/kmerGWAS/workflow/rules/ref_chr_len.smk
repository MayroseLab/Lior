rule ref_chr_len:
    """
    Create a file with reference
    genome chromosome lengths
    """
    input:
        os.path.join(out_dir, 'genomes', '%s.fasta' % ref_genome)
    output:
        os.path.join(out_dir, 'genomes', '%s.chr.len' % ref_genome)
    log:
        os.path.join(logs_dir, 'ref_chr_len', 'ref_chr_len.log')
    conda:
        os.path.join(envs_dir, 'bioawk.yaml')
    shell:
        """
        bioawk -c fastx '{{ print $name, length($seq) }}' {input} > {output}
        """
