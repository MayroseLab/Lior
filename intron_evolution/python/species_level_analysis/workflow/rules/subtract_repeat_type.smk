rule subtract_repeat_type:
    """
    Subtract each of the repeat
    types from "all repeats".
    This is helpful for assessing
    the effect of each repeat type
    Also subtract repet type from
    introns and PIRs.
    """
    input:
        all_repeats_bed = os.path.join(out_dir, 'per_species', '{species}', 'repeats', 'all.merge.bed'),
        rep_type_bed = os.path.join(out_dir, 'per_species', '{species}', 'repeats', '{rep_type}.merge.bed'),
        introns_bed = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.bed'),
        introns_flank_size_limit_bed = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.flank_size_limit.bed'),
        intergenic_bed = os.path.join(out_dir, 'per_species', '{species}', 'PIR.bed')
    output:
        repeats = os.path.join(out_dir, 'per_species', '{species}', 'repeats', 'all-{rep_type}.merge.bed'),
        introns = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns-{rep_type}.bed'),
        introns_flank_size_limit = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.flank_size_limit-{rep_type}.bed'),
        intergenic = os.path.join(out_dir, 'per_species', '{species}', 'PIR-{rep_type}.bed')
    log:
        os.path.join(logs_dir, 'subtract_repeat_type', '{species}.{rep_type}.subtract_repeat_type.log')
    conda:
        os.path.join(envs_dir, 'bedtools.yml')
    resources:
        mem_mb = 4000
    shell:
        """
        if [ {wildcards.rep_type} == "all" ]
        then
            touch {output} 
        else
            bedtools subtract -a {input.all_repeats_bed} -b {input.rep_type_bed} > {output.repeats}
            bedtools subtract -a {input.introns_bed} -b {input.rep_type_bed} > {output.introns}
            bedtools subtract -a {input.introns_flank_size_limit_bed} -b {input.rep_type_bed} > {output.introns_flank_size_limit}
            bedtools subtract -a {input.intergenic_bed} -b {input.rep_type_bed} > {output.intergenic}
        fi
        """
