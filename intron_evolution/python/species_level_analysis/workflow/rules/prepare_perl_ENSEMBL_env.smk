rule prepare_perl_ENSEMBL_env:
    """
    Manually prepare required libraries
    for using the ENSEMBL Perl API.
    (cannot be done with conda)
    """
    output:
        os.path.join(out_dir, 'perl_ensembl_env.done')
    log:
        os.path.join(logs_dir, 'prepare_perl_ENSEMBL_env', 'prepare_perl_ENSEMBL_env.log')
    conda:
        os.path.join(envs_dir, 'bioperl.yml')
    shell:
        """
        cd $CONDA_PREFIX
        mkdir src
        cd src
        wget https://ftp.ensembl.org/pub/ensembl-api.tar.gz
        wget https://github.com/bioperl/bioperl-live/archive/release-1-6-924.zip
        tar zxvf ensembl-api.tar.gz
        unzip release-1-6-924.zip
        mv bioperl-live-release-1-6-924 bioperl-1.6.924
        touch {output}
        """
