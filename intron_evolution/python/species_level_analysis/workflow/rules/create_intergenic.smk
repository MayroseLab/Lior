rule create_intergenic:
    """
    Create a BED file of intergenic
    regions - those flanking genes
    but not overlapping any other gene
    """
    input:
        gff = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.gff3'),
        genome_file = os.path.join(out_dir, 'per_species', '{species}', 'genome.sm.genome')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'intergenic.bed')
    log:
        os.path.join(logs_dir, 'create_intergenic', '{species}.create_intergenic.log')
    conda:
        os.path.join(envs_dir, 'bedtools.yml')
    params:
        flank_ize = config['flank_size']
    resources:
        mem_mb = 4000
    shell:
        """
        # create genes bed
        awk '$3 == "gene" {{print $1"\t"$4-1"\t"$5}}' {input.gff} > {input.gff}.genes.bed
        # create flank bed
        bedtools flank -i {input.gff}.genes.bed -g {input.genome_file} -b {params.flank_ize} | sort -k1,1 -k2,2n > {input.gff}.genes.flank.bed
        # merge overlapping regions
        bedtools merge -i {input.gff}.genes.flank.bed > {input.gff}.genes.flank.merge.bed
        # subtract genes
        bedtools subtract -a {input.gff}.genes.flank.merge.bed -b {input.gff}.genes.bed > {output}
        """
