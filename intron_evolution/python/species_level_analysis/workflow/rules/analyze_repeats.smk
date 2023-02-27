rule analyze_repeats:
    """
    Analyze and extract stats
    about genomic repeats
    """
    input:
        introns_bed = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.bed'),
        introns_flank_size_limit_bed = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.flank_size_limit.bed'),
        intergenic_bed = os.path.join(out_dir, 'per_species', '{species}', 'PIR.bed'),
        repeats_bed = os.path.join(out_dir, 'per_species', '{species}', 'repeats', '{rep_type}.merge.bed'),
        genome_file = os.path.join(out_dir, 'per_species', '{species}', 'genome.sm.genome')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'repeats.{rep_type}.stats')
    log:
        os.path.join(logs_dir, 'analyze_repeats', '{species}.{rep_type}.analyze_repeats.log')
    conda:
        os.path.join(envs_dir, 'bedtools.yml')
    resources:
        mem_mb = 6000
    shell:
        """
        tot_repeats=$(awk '{{SUM+=$3-$2}}END{{print SUM}}' {input.repeats_bed})
        tot_int=$(awk '{{SUM+=$3-$2}}END{{print SUM}}' {input.introns_bed})
        tot_asm=$(awk '{{SUM+=$2}}END{{print SUM}}' {input.genome_file})
        tot_intergenic=$(awk '{{SUM+=$3-$2}}END{{print SUM}}' {input.intergenic_bed})

        perc_rep_genome=$(echo "$tot_repeats/$tot_asm*100" | bc -l)
        rep_in_introns=$(bedtools intersect -a {input.introns_bed} -b {input.repeats_bed} | awk '{{SUM+=$3-$2}}END{{print SUM}}')
        perc_rep_introns=$(echo "$rep_in_introns/$tot_int*100" | bc -l)
        rep_in_inter=$(bedtools intersect -a {input.intergenic_bed} -b {input.repeats_bed} | awk '{{SUM+=$3-$2}}END{{print SUM}}')
        perc_rep_intergenic=$(echo "$rep_in_inter/$tot_intergenic*100" | bc -l)

        mean_rep_len_introns=$(bedtools intersect -a {input.repeats_bed} -b {input.introns_bed} -f 1 -wa | awk '{{SUM+=$3-$2; N+=1}}END{{ if (N==0) {{print 0}} else {{print SUM/N}}}}')
        mean_rep_len_introns_limit=$(bedtools intersect -a {input.repeats_bed} -b {input.introns_flank_size_limit_bed} -f 1 -wa | awk '{{SUM+=$3-$2; N+=1}}END{{ if (N==0) {{print 0}} else {{print SUM/N}}}}')
        mean_rep_len_intergenic=$(bedtools intersect -a {input.repeats_bed} -b {input.intergenic_bed} -f 1 -wa | awk '{{SUM+=$3-$2; N+=1}}END{{ if (N==0) {{print 0}} else {{print SUM/N}}}}')

        echo -e "{wildcards.species}\t{wildcards.rep_type}\t$tot_asm\t$tot_repeats\t$tot_int\t$tot_intergenic\t$perc_rep_genome\t$perc_rep_introns\t$perc_rep_intergenic\t$mean_rep_len_introns\t$mean_rep_len_introns_limit\t$mean_rep_len_intergenic" > {output}
        """
