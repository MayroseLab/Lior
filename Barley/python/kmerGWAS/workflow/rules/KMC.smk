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

rule count_canonized_kmers:
    """
    Count k-mers using canonization
    """
    input:
        files_list = os.path.join(out_dir, 'per_sample', '{sample}', 'KMC.lst'),
    output:
        suf = temp(os.path.join(out_dir, 'per_sample', '{sample}','{sample}_canon.kmc_suf')),
        pre = temp(os.path.join(out_dir, 'per_sample', '{sample}','{sample}_canon.kmc_pre'))
    params:
        k = config['k'],
        min_count = config['min_count'],
        file_type = get_file_type,
        out_dir = os.path.join(out_dir, 'per_sample', '{sample}'),
    threads:
        config['cpu']
    resources:
        mem_gb = config['mem_gb']
    conda:
        os.path.join(envs_dir, 'KMC.yaml')
    log:
        os.path.join(logs_dir, 'KMC', '{sample}.count_canonized_kmers.log')
    shell:
        """
        cd {params.out_dir}
        kmc -k{params.k} -m24 -ci{params.min_count} -t{threads} -{params.file_type} @{input.files_list} {wildcards.sample}_canon {params.out_dir} --opt-out-size
        """

rule count_uncanonized_kmers:
    """
    Count k-mers without canonization
    """
    input:
        files_list = os.path.join(out_dir, 'per_sample', '{sample}', 'KMC.lst'),
    output:
        suf = temp(os.path.join(out_dir, 'per_sample', '{sample}','{sample}_all.kmc_suf')),
        pre = temp(os.path.join(out_dir, 'per_sample', '{sample}','{sample}_all.kmc_pre'))
    params:
        k = config['k'],
        file_type = get_file_type,
        out_dir = os.path.join(out_dir, 'per_sample', '{sample}')
    threads:
        config['cpu']
    resources:
        mem_gb = config['mem_gb']
    conda:
        os.path.join(envs_dir, 'KMC.yaml')
    log:
        os.path.join(logs_dir, 'KMC', '{sample}.count_uncanonized_kmers.log')
    shell:
        """
        cd {params.out_dir}
        kmc -k{params.k} -m24 -ci0 -b -t{threads} -{params.file_type} @{input.files_list} {wildcards.sample}_all {params.out_dir} --opt-out-size
        """

rule combine_kmer_counts:
    """
    Combine canonized and uncanonized k-mer counts
    and add strand information
    """
    input:
        os.path.join(out_dir, 'per_sample', '{sample}','{sample}_canon.kmc_suf'),
        os.path.join(out_dir, 'per_sample', '{sample}','{sample}_canon.kmc_pre'),
        os.path.join(out_dir, 'per_sample', '{sample}','{sample}_all.kmc_suf'),
        os.path.join(out_dir, 'per_sample', '{sample}','{sample}_all.kmc_pre')
    output:
        os.path.join(out_dir, 'per_sample', '{sample}','{sample}_kmers_with_strand')
    params:
        k = config['k'],
        add_strand_script = os.path.join(config['voichek_code_dir'], 'bin/kmers_add_strand_information'),
        out_dir = os.path.join(out_dir, 'per_sample', '{sample}')
    resources:
        mem_gb = config['mem_gb']
    conda:
        os.path.join(envs_dir, 'gcc.yaml')
    log:
        os.path.join(logs_dir, 'KMC', '{sample}.combine_kmer_counts.log')
    shell:
        """
        cd {params.out_dir}
        {params.add_strand_script} -c {wildcards.sample}_canon -n {wildcards.sample}_all -k {params.k} -o {output}
        """
