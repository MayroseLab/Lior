ensembl_vert_ftp = "ftp.ensembl.org"
ensembl_genomes_ftp = "ftp.ensemblgenomes.ebi.ac.uk"

rule download_genome:
    """
    Download genome sequence
    (soft-masked) in FASTA format
    """
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'genome.sm.fasta')
    log:
        os.path.join(logs_dir, 'download_genome', '{species}.download_genome.log')
    run:
        species = wildcards.species
        group = species_df.loc[species]['group']
        if group == 'vertebrates':
            ftp = ftputil.FTPHost(ensembl_vert_ftp, "anonymous", "")
            ftp.chdir(f'pub/current_fasta/{species}/dna')           
        else:
            ftp = ftputil.FTPHost(ensembl_genomes_ftp, "anonymous", "")
            ftp.chdir(f'pub/current/{group}/fasta/{species}/dna')
        all_files = ftp.listdir(ftp.curdir)
        # Check if "primary assembly" exists
        prim_asm = [f for f in all_files if f.endswith('_sm.primary_assembly.fa.gz')]
        if len(prim_asm) > 0:
            asm_files = prim_asm
        # if not, take "toplevel"
        else:
            asm_files = [f for f in all_files if f.endswith('.dna_sm.toplevel.fa.gz')]
        if len(asm_files) != 1:
            print(asm_files)
            sys.exit(f"Can't determine soft-masked genome for species {species}")
        sm_genome = asm_files[0]
        ftp.download(sm_genome, output[0] + '.gz')
        os.system('gzip -d ' + output[0] + '.gz')
