
rule remove_alt_splicing:
    """
    Remove alternative splice
    variants for all genes, leaving
    only the canonical transcript.
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.gff3')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.gff3')
    log:
        os.path.join(logs_dir, 'download_annotation', '{species}.remove_alt_splicing.log')
    conda:
        os.path.join(envs_dir, 'gffutils.yml')
    params:
        remove_alt_splicing_script = os.path.join(scripts_dir, 'remove_alt_splicing_from_gff.py')
    shell:
        """
        python {params.remove_alt_splicing_script} {input} {output}
        """
