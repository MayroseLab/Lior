from snakemakeUtils import *

def init(): 
    #load_info_file
    config['samples_info'] = SampleInfoReader.sample_table_reader(filename=config['samples_info_file'], 
                delimiter='\t', key_name='sample', col_names=['path'])

init()

onstart:
    write_config_file(config)

onsuccess:
    print("%s pipeline finished, no error" % config['name'])
    shell("cat {log} >> %s/run_log.txt" % config["out_dir"])

onerror:
    print("%s pipeline failed" % config['name'])
    shell("cat {log} >> %s/run_log.txt" % config["out_dir"])

#------------------------------------
#                RULES              |
#------------------------------------

localrules: all, prep_maker_liftover_configs, prep_maker_annotation_configs 

rule all:
    input:
        expand(config["out_dir"] + "/per_sample/{sample}/Annotation_QA/{sample}.QA_report.tsv", sample=config['samples_info']) 
	
def get_sample(wildcards):
    return config['samples_info'][wildcards.sample]['path']

rule prep_maker_liftover_configs:
    input:
        get_sample
    params:
        templates=config["config_templates"],
        official_transcripts=config["official_transcripts_fasta"]
    output:
        bopts=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_bopts.ctl",
        opts=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_opts.ctl",
        exe=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_exe.ctl"
    shell:
        """
        cp {params.templates}/maker_bopts.ctl {output.bopts}
        cp {params.templates}/maker_exe.ctl {output.exe}
        sed -e 's|<genome_fasta>|{input}|' -e 's|<transcripts_fasta>|{params.official_transcripts}|' {params.templates}/maker_opts_liftover.ctl > {output.opts}
        """

rule maker_liftover:
    input:
        bopts=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_bopts.ctl",
        opts=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_opts.ctl",
        exe=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_exe.ctl"
    output:
        config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/liftover.done",
        config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/{sample}_genome_liftover.maker.output/{sample}_genome_liftover_master_datastore_index.log"
    params:
        run_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover",
        base_name="{sample}_genome_liftover",
        nodes=config['nodes'],
        ppn=config['ppn']
    threads:
        config['ppn'] * config['nodes']
    conda:
        "conda_env/maker-env.yaml"
    shell:
        """
        cd {params.run_dir}
        module load miniconda/miniconda2-4.5.4-MakerMPI
        mpirun -n {threads} maker -base {params.base_name}
        if [ -f {params.base_name}.maker.output/{params.base_name}_master_datastore_index.log ] && [ `grep STARTED {params.base_name}.maker.output/{params.base_name}_master_datastore_index.log | wc -l` == `grep FINISHED {params.base_name}.maker.output/{params.base_name}_master_datastore_index.log | wc -l` ]; then touch {output}; fi
        """

rule create_liftover_gff:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/{sample}_genome_liftover.maker.output/{sample}_genome_liftover_master_datastore_index.log"
    output:
        config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/{sample}_genome_liftover.maker.output/{sample}_genome_liftover.all.gff"
    threads:
        1
    params:
        nodes=1,
        ppn=1,
        liftover_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/{sample}_genome_liftover.maker.output"
    conda:
        "conda_env/maker-env.yaml"
    shell:
        """
        cd {params.liftover_dir}
        gff3_merge -d {input} -g -n
        """

rule prep_maker_annotation_configs:
    input:
        genome=get_sample,
        liftover_gff=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/{sample}_genome_liftover.maker.output/{sample}_genome_liftover.all.gff"
    params:
        templates=config["config_templates"],
        transcripts=config['annotation_transcripts_fasta'],
        proteins=config['annotation_proteins_fasta'],
        external_gff=config['external_gff']
    output:
        bopts=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_bopts.ctl",
        opts=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_opts.ctl",
        exe=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_exe.ctl"
    run:
        shell("cp {params.templates}/maker_bopts.ctl {output.bopts}")
        shell("cp {params.templates}/maker_exe.ctl {output.exe}")
        if params.external_gff:
            shell("sed -e 's|<genome_fasta>|{input.genome}|' -e 's|<transcripts_fasta>|{params.transcripts}|' -e 's|<proteins_fasta>|{params.proteins}|' -e 's|<gene_models_gff>|{input.liftover_gff},{params.external_gff}|' {params.templates}/maker_opts.ctl > {output.opts}")
        else:
            shell("sed -e 's|<genome_fasta>|{input.genome}|' -e 's|<transcripts_fasta>|{params.transcripts}|' -e 's|<proteins_fasta>|{params.proteins}|' -e 's|<gene_models_gff>|{input.liftover_gff}|' {params.templates}/maker_opts.ctl > {output.opts}")

rule maker_annotation:
    input:
        bopts=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_bopts.ctl",
        opts=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_opts.ctl",
        exe=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_exe.ctl"
    output:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/annotation.done",
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}_master_datastore_index.log"
    params:
        run_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation",
        base_name="{sample}",
        nodes=config['nodes'],
        ppn=config['ppn']
    threads:
        config['ppn'] * config['nodes']
    conda:
        "conda_env/maker-env.yaml"
    shell:
        """
        cd {params.run_dir}
        module load miniconda/miniconda2-4.5.4-MakerMPI
        mpirun -n {threads} maker -base {params.base_name}
        if [ -f {params.base_name}.maker.output/{params.base_name}_master_datastore_index.log ] && [ `grep STARTED {params.base_name}.maker.output/{params.base_name}_master_datastore_index.log | wc -l` == `grep FINISHED {params.base_name}.maker.output/{params.base_name}_master_datastore_index.log | wc -l` ]; then touch {output}; fi
        """
rule create_annotation_gff:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}_master_datastore_index.log"
    output:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}.all.gff"
    threads:
        1
    params:
        nodes=1,
        ppn=1,
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output"
    conda:
        "conda_env/maker-env.yaml"
    shell:
        """
        cd {params.annotation_dir}
        gff3_merge -d {input} -n
        """

rule create_annotation_fasta:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}_master_datastore_index.log"
    output:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}.all.maker.proteins.fasta"
    threads:
        1
    params:
        nodes=1,
        ppn=1,
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output"
    conda:
        "conda_env/maker-env.yaml"
    shell:
        """
        cd {params.annotation_dir}
        fasta_merge -d {input}
        """

rule proteins_busco:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}.all.maker.proteins.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/run_BUSCO/full_table_BUSCO.tsv"
    params:
        busco_exe=config['busco_exe'],
        busco_set=config['busco_set'],
        augustus_conf=config['augustus_conf'],
        augustus_bin=config['augustus_bin'],
        out_dir=config["out_dir"] + "/per_sample/{sample}/",
        nodes=1,
        ppn=config['ppn']
    threads:
        config['ppn']
    shell:
        """
        export AUGUSTUS_CONFIG_PATH={params.augustus_conf}
        export PATH={params.augustus_bin}:$PATH
        cd {params.out_dir}
        python {params.busco_exe} --in {input} --out BUSCO --lineage_path {params.busco_set} --cpu {threads} --mode proteins -f
        """

rule proteins_interproscan:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}.all.maker.proteins.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/interProScan/{sample}.interProScan.tsv"
    params:       
        nodes=1,
        ppn=1,
        ips_dir=config["out_dir"] + "/per_sample/{sample}/interProScan/"
    shell:
        """
        module load interproscan/5.32-71
        module load java/jdk1.8.25
        interproscan.sh -i {input} -b {params.ips_dir}/{wildcards.sample}.interProScan -t p -dp -pa -appl Pfam,ProDom,SuperFamily,PIRSF --goterms --iprlookup -f tsv,gff3 -T {params.ips_dir}
        """

rule proteins_blast:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}.all.maker.proteins.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/Blast/{sample}.blast.tsv"
    conda:
        "conda_env/blast.yaml"
    params:
        nodes=1,
        ppn=config['ppn'],
        blast_qa_db=config['blast_qa_db']
    threads:
        config['ppn']
    shell:
        """
        blastp -query {input} -db {params.blast_qa_db} -evalue 1e-5 -outfmt \"6 qseqid sseqid pident qlen length mismatch gapopen evalue bitscore qcovs\" -num_threads {threads} -max_target_seqs 1 -out {output}
        """

rule annotation_qa:
    input:
        gff=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}.all.gff",
        busco=config["out_dir"] + "/per_sample/{sample}/run_BUSCO/full_table_BUSCO.tsv",
        ips=config["out_dir"] + "/per_sample/{sample}/interProScan/{sample}.interProScan.tsv",
        blast=config["out_dir"] + "/per_sample/{sample}/Blast/{sample}.blast.tsv"
    output:
        config["out_dir"] + "/per_sample/{sample}/Annotation_QA/{sample}.QA_report.tsv"
    params:
        nodes=1,
        ppn=1,
        qa_script=config['qa_script']
    conda:
        "conda_env/annotation_qa.yaml"
    shell:
        """
        python {params.qa_script} {input.gff} {output} --busco_result {input.busco} --blast_result {input.blast} --ips_result {input.ips}
        """
