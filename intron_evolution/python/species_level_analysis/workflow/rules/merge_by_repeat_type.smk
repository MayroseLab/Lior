rule merge_by_repeat_type:
    """
    Takes a BED file with repeat
    features, splits into multiple
    BEDs based on repeat type
    (column 6) and bedtools-merges
    features in each file. Also merges
    all overlapping repeats regardless
    of type.
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'repeats', 'repeats.bed')
    output:
        merged_type = os.path.join(out_dir, 'per_species', '{species}', 'repeats', '{rep_type}.merge.bed'),
        all_minus_type = os.path.join(out_dir, 'per_species', '{species}', 'repeats', 'all-{rep_type}.merge.bed')
    log:
        os.path.join(logs_dir, 'merge_by_repeat_type', '{species}.{rep_type}.merge_by_repeat_type.log')
    conda:
        os.path.join(envs_dir, 'bedtools.yml')
    resources:
        mem_mb = 4000
    shell:
        """
        if [ {wildcards.rep_type} == "all" ]
        then
            cut -f1,2,3 {input} | sort -k1,1 -k2,2n | bedtools merge > {output.merged_type}
            touch {output.all_minus_type}
        else
            awk '$6 == "{wildcards.rep_type}"' {input} | cut -f1,2,3 | sort -k1,1 -k2,2n | bedtools merge > {output.merged_type}
            awk '$6 != "{wildcards.rep_type}"' {input} | cut -f1,2,3 | sort -k1,1 -k2,2n | bedtools merge > {output.all_minus_type}
        fi
        """
