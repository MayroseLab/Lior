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
        os.path.join(out_dir, 'per_species', '{species}', 'repeats', '{rep_type}.merge.bed')
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
            sort -k1,1 -k2,2n {input} | bedtools merge > {output}
        else
            awk '$6 == "{wildcards.rep_type}"' {input} | sort -k1,1 -k2,2n | bedtools merge > {output}
        fi
        """
