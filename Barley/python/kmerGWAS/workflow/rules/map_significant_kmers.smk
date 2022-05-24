rule map_kmers:
    """
    Map significant k-mer sequences
    to a genome using the BWA-Backtrack
    algorithm
    """
    input:
        kmers_fasta = os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per.fasta'),
        genome_fasta = os.path.join(out_dir, 'genomes', '{genome}.fasta'),
        genome_index = os.path.join(out_dir, 'genomes', '{genome}.fasta.bwt')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_vs_{genome}.sai')
    log:
        os.path.join(logs_dir, 'map_significant_kmers', '{phenotype}.{genome}.map_kmers.log')
    conda:
        os.path.join(envs_dir, 'bwa.yaml')
    resources:
        mem_gb = 16
    shell:
        """
        bwa aln {input.genome_fasta} {input.kmers_fasta} > {output}
        """

rule sai_to_sam:
    """
    Convert BWA-Backtrack binary
    output to SAM
    """
    input:
        sai = os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_vs_{genome}.sai'),
        kmers_fasta = os.path.join(out_dir, 'all_samples', '{phenotype}', 'pass_threshold_5per.fasta'),
        genome_fasta = os.path.join(out_dir, 'genomes', '{genome}.fasta')
    output:
        os.path.join(out_dir, 'all_samples', '{phenotype}', 'map_signif_kmers', 'pass_threshold_5per_vs_{genome}.sam')
    log:
        os.path.join(logs_dir, 'map_significant_kmers', '{phenotype}.{genome}.sai_to_sam.log')
    conda:
        os.path.join(envs_dir, 'bwa.yaml')
    resources:
        mem_gb = 32
    shell:
        """
        bwa samse {input.genome_fasta} {input.sai} {input.kmers_fasta} -f {output}
        """
