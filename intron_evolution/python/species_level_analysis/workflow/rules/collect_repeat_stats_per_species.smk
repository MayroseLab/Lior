rule collect_repeat_stats_per_species:
    """
    Collect repeat stats for
    a given species for different
    repeat types
    """
    input:
        expand(os.path.join(out_dir, 'per_species', '{{species}}', 'repeats.{rep_type}.stats'), rep_type=['all'] + repeat_types + ['all-%s' % rt for rt in repeat_types])
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'repeats.stats')
    log:
        os.path.join(logs_dir, 'collect_repeat_stats_per_species', '{species}.collect_repeat_stats_per_species.log')
    shell:
        """
        cat {input} > {output}
        """
