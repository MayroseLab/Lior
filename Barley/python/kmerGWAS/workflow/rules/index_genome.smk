rule copy_genome:
    """
    Copy all genomes fasta
    to one place
    """
    input:
        get_genome
    output:
        os.path.join(out_dir, 'genomes', '{genome}.fasta')
    log:
        os.path.join(logs_dir, 'index_genome', '{genome}.copy_genome.log')
    shell:
        """
        cp {input} {output}
        """

rule index_genome:
    """
    BWA index of all high uqality
    and low quality genome assemblies
    """
    input:
        os.path.join(out_dir, 'genomes', '{genome}.fasta')
    output:
        multiext(os.path.join(out_dir, 'genomes', '{genome}.fasta'), '.bwt', '.pac', '.ann', '.amb', '.sa')
    log:
        os.path.join(logs_dir, 'index_genome', '{genome}.index_genome.log')
    conda:
        os.path.join(envs_dir, 'bwa.yaml')
    resources:
        mem_gb = 16
    shell:
        """
        bwa index {input}
        """
