import os

rule make_files_list:
    """
    Make files list to be used by KMC
    """
    input:
        get_reads
    output:
        os.path.join(out_dir, 'per_sample', '{sample}', 'KMC.lst')
    log:
        os.path.join(logs_dir, 'KMC', '{sample}.make_files_list.log')
    shell:
        """
        echo {input} | tr ' ' '\n' > {output}
        """

rule run_KMC:
    """
    Count k-mers and create histogram
    """
    input:
        files_list = os.path.join(out_dir, 'per_sample', '{sample}', 'KMC.lst'),
    output:
        suf = temp(os.path.join(out_dir, 'per_sample', '{sample}','{sample}.kmc_suf')),
        pre = temp(os.path.join(out_dir, 'per_sample', '{sample}','{sample}.kmc_pre')),
        suf_sort = os.path.join(out_dir, 'per_sample', '{sample}','{sample}_sort.kmc_suf'),
        pre_sort = os.path.join(out_dir, 'per_sample', '{sample}','{sample}_sort.kmc_pre')
    params:
        k = config['k'],
        min_count = config['min_count'],
        max_count = config['max_count'],
        file_type = get_file_type,
        out_dir = os.path.join(out_dir, 'per_sample', '{sample}'),
    threads:
        config['cpu']
    resources:
        mem_gb = config['mem_gb']
    conda:
        os.path.join(envs_dir, 'KMC.yaml')
    log:
        os.path.join(logs_dir, 'KMC', '{sample}.run_KMC.log')
    shell:
        """
        cd {params.out_dir}
        kmc -k{params.k} -m24 -ci{params.min_count} -cx{params.max_count} -cs{params.max_count} -t{threads} -{params.file_type} @{input.files_list} {wildcards.sample} {params.out_dir} --opt-out-size
        kmc_tools transform {wildcards.sample} sort {wildcards.sample}_sort
        """

rule create_matrixer_list:
    input:
        pre = expand(os.path.join(out_dir, 'per_sample', '{sample}','{sample}_sort.kmc_pre'), sample=samples),
        suf = expand(os.path.join(out_dir, 'per_sample', '{sample}','{sample}_sort.kmc_suf'), sample=samples)
    output:
        os.path.join(out_dir, 'all_samples', 'kmer_matrix.lst')
    log:
        os.path.join(logs_dir, 'KMC', 'all_samples.create_matrixer_list.log')
    shell:
        """
        echo "{input.pre}" | sed 's/\.kmc_pre//g' | tr ' ' '\n' > {output}
        """

rule create_kmer_matrix:
    input:
        os.path.join(out_dir, 'all_samples', 'kmer_matrix.lst')
    output:
        os.path.join(out_dir, 'all_samples', 'kmer_matrix.tsv')
    params:
        matrixer_exec = config['matrixer_exec']
    resources:
        mem_gb = config['mem_gb']
    conda:
        os.path.join(envs_dir, 'gcc.yaml')
    log:
        os.path.join(logs_dir, 'KMC', 'all_samples.create_kmer_matrix.log')
    shell:
        """
        {params.matrixer_exec} {input} > {output}
        """
    
