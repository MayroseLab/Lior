rule calculate_distance_matrix:
    """
    """
    input:
        expand(os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.gff3'), species=species_list)
    output:
        os.path.join(out_dir, 'all_species', 'distance_matrix.tsv')
    log:
        os.path.join(logs_dir, 'calculate_distance_matrix', 'all_species.calculate_distance_matrix.log')
    params:
        matrix_script = os.path.join(scripts_dir, 'calculate_distance_matrix.py')
    conda:
        os.path.join(envs_dir, 'python_analyze.yml')
    shell:
        """
        python {params.matrix_script} {output} {input}
        """
