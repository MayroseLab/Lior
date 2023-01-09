def get_group(wc):
    species = wc.species
    return species_df.loc[species]['group']

rule download_repeats:
    """
    Download repeats annotation
    from ENSEMBL API in BED format
    """
    input:
        os.path.join(out_dir, 'perl_ensembl_env.done')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'repeats', 'repeats.bed')
    log:
        os.path.join(logs_dir, 'download_repeats', '{species}.download_repeats.log')
    conda:
        os.path.join(envs_dir, 'bioperl.yml')
    params:
        download_script = os.path.join(scripts_dir, 'download_repeats_from_ENSEMBL_API.pl'),
        group = get_group
    shell:
        """
        PERL5LIB=${{CONDA_PREFIX}}/src/bioperl-1.6.924
        PERL5LIB=${{PERL5LIB}}:${{CONDA_PREFIX}}/src/ensembl/modules
        PERL5LIB=${{PERL5LIB}}:${{CONDA_PREFIX}}/src/ensembl-compara/modules
        PERL5LIB=${{PERL5LIB}}:${{CONDA_PREFIX}}/src/ensembl-variation/modules
        PERL5LIB=${{PERL5LIB}}:${{CONDA_PREFIX}}/src/ensembl-funcgen/modules
        PERL5LIB=${{PERL5LIB}}:${{CONDA_PREFIX}}/src/ensembl-metadata/modules
        export PERL5LIB
        perl {params.download_script} {wildcards.species} {params.group} {output}
        """
