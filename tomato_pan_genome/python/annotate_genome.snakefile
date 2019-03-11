from snakemakeUtils import *
from time import time

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

localrules: all, prep_maker_liftover_configs, prep_maker_annotation_configs, maker_liftover_chunks, maker_annotation_chunks, prep_chunks_tsv, prep_liftover_yaml, prep_annotation_yaml

rule all:
    input:
        expand(config["out_dir"] + "/per_sample/{sample}/Annotation_QA/{sample}.QA_report.tsv", sample=config['samples_info']) 

def get_sample(wildcards):
    return config['samples_info'][wildcards.sample]['path']

if config['maker_parallel'] == "chunks":
    
    rule break_genome_to_chunks:
        input:
            get_sample
        output:
            config["out_dir"] + "/per_sample/{sample}/genome_chunks/chunks.done"
        params:
            n_chunks=config['n_chunks'],
            break_script=config['break_to_chunks_script'],
            chunks_out_dir=config["out_dir"] + "/per_sample/{sample}/genome_chunks",
            nodes=1,
            ppn=1
        conda:
            "conda_env/biopython.yaml"
        shell:
            """
            #mkdir {params.chunks_out_dir}
            python {params.break_script} {input} {params.chunks_out_dir} --n_chunks {params.n_chunks} --break_mode
            if (( `ls -1 {params.chunks_out_dir} | wc -l` >= {params.n_chunks} )); then touch {output}; fi
            """

    rule prep_chunks_tsv:
        input:
            config["out_dir"] + "/per_sample/{sample}/genome_chunks/chunks.done"
        output:
            config["out_dir"] + "/per_sample/{sample}/genome_chunks/chunks.tsv"
        params:
            chunks_dir=config["out_dir"] + "/per_sample/{sample}/genome_chunks"
        shell:
            """
            echo -e "chunk\tpath" > {output}
            for f in `ls -1 {params.chunks_dir}/*.fasta`; do chunkName=`basename $f | sed 's/\.fasta//'`; echo -e "$chunkName\t$f" >> {output}; done
            """

    rule prep_liftover_yaml:
        input:
            config["out_dir"] + "/per_sample/{sample}/genome_chunks/chunks.tsv"
        output:
            config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/config.yaml"
        params:
            liftover_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover",
            templates_dir=config["liftover_config_templates"],
            coord_conversion_script=config['coord_conversion_script'],
            config_edit_script=config['config_edit_script'],
            official_transcripts_fasta=config['official_transcripts_fasta']
        shell:
            """               
            echo "name: MAKER_wrapper" >> {output}
            echo "chunks_info_file: {input}" >> {output}
            echo "out_dir: {params.liftover_dir}" >> {output}
            echo "config_templates: {params.templates_dir}" >> {output}
            echo "coord_conversion_script: {params.coord_conversion_script}" >> {output}
            echo "config_edit_script: {params.config_edit_script}" >> {output}
            echo config_kv_pairs: est={params.official_transcripts_fasta} >> {output}
            """

    rule maker_liftover_chunks:
        input:
            config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/config.yaml"
        output:
            config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker.genes.convert.gff"
        params:
            run_maker_in_chunks_snakefile=config['run_maker_in_chunks_snakefile'],
            queue=config['queue'],
            jobs=config['jobs']//len(config['samples_info']),
            liftover_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover"
        shell:
            """
            cd {params.liftover_dir}
            snakemake -s {params.run_maker_in_chunks_snakefile} --configfile {input} --cluster "qsub -N {wildcards.sample}_chunk_liftover -q {params.queue}" -j {params.jobs} --latency-wait 60 --restart-times 3
            """

elif config['maker_parallel'] == "mpi":

    rule prep_maker_liftover_configs:
        input:
            get_sample
        params:
            templates=config["config_templates"] + "/liftover",
            official_transcripts=config["official_transcripts_fasta"]
        output:
            bopts=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_bopts.ctl",
            opts=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_opts.ctl",
            exe=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_exe.ctl"
        log:
            stdout=config["out_dir"] + "/logs/{sample}_prep_maker_liftover_configs_" + str(time()) + ".out",
            stderr=config["out_dir"] + "/logs/{sample}_prep_maker_liftover_configs_" + str(time()) + ".err"
        shell:
            """
            cp {params.templates}/maker_bopts.ctl {output.bopts} > {log.stdout} 2> {log.stderr}
            cp {params.templates}/maker_exe.ctl {output.exe} >> {log.stdout} 2>> {log.stderr}
            sed -e 's|<genome_fasta>|{input}|' -e 's|<transcripts_fasta>|{params.official_transcripts}|' {params.templates}/maker_opts.ctl > {output.opts} 2>> {log.stderr}
            """
    
    rule maker_liftover:
        input:
            bopts=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_bopts.ctl",
            opts=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_opts.ctl",
            exe=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker_exe.ctl"
        output:
            config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/liftover.done"
        log:
            index=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/{sample}_genome_liftover.maker.output/{sample}_genome_liftover_master_datastore_index.log",
            stdout=config["out_dir"] + "/logs/{sample}_maker_liftover_" + str(time()) + ".out",
            stderr=config["out_dir"] + "/logs/{sample}_maker_liftover_" + str(time()) + ".err"
        params:
            run_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover",
            base_name="{sample}_genome_liftover",
            nodes=config['nodes'],
            ppn=config['ppn']
        threads:
            config['ppn'] * config['nodes']
    #    conda:
    #        "conda_env/maker-env.yaml"
        shell:
            """
            cd {params.run_dir}
            module load miniconda/miniconda2-4.5.4-MakerMPI
            mpirun -n {threads} -env I_MPI_FABRICS tcp maker -base {params.base_name} > {log.stdout} 2> {log.stderr}
            if [ -f {log.index} ] && [ `grep STARTED {log.index} | wc -l` == `grep FINISHED {log.index} | wc -l` ]; then touch {output}; fi
            """
    
    rule create_liftover_gff:
        input:
            config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/{sample}_genome_liftover.maker.output/{sample}_genome_liftover_master_datastore_index.log"
        output:
            config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/{sample}_genome_liftover.maker.output/{sample}_genome_liftover.all.gff"
        log:
            stdout=config["out_dir"] + "/logs/{sample}_create_liftover_gff_" + str(time()) + ".out",
            stderr=config["out_dir"] + "/logs/{sample}_create_liftover_gff_" + str(time()) + ".err"
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
            gff3_merge -d {input} -g -n > {log.stdout} 2> {log.stderr}
            """

if config['maker_parallel'] == "chunks":

    rule prep_annotation_yaml:
        input:
            chunks_tsv=config["out_dir"] + "/per_sample/{sample}/genome_chunks/chunks.tsv",
            liftover_gff=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/maker.genes.convert.gff"
        output:
            config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/config.yaml"
        params:
            annotation_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation",
            templates_dir=config["annotation_config_templates"],
            coord_conversion_script=config['coord_conversion_script'],
            config_edit_script=config['config_edit_script'],
            transcripts=config['annotation_transcripts_fasta'],
            proteins=config['annotation_proteins_fasta'],
            external_gff=config['external_gff']
        shell:
            """               
            echo "name: MAKER_wrapper" >> {output}
            echo "chunks_info_file: {input.chunks_tsv}" >> {output}
            echo "out_dir: {params.annotation_dir}" >> {output}
            echo "config_templates: {params.templates_dir}" >> {output}
            echo "coord_conversion_script: {params.coord_conversion_script}" >> {output}
            echo "config_edit_script: {params.config_edit_script}" >> {output}
            echo config_kv_pairs: est={params.transcripts} protein={params.proteins} pred_gff={input.liftover_gff} >> {output}
            """

    rule maker_annotation_chunks:
        input:
            config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/config.yaml"
        output:
            config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.genes.convert.gff",
            config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.all.convert.gff",
            config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.proteins.fasta",
            config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.transcripts.fasta"
        params:
            run_maker_in_chunks_snakefile=config['run_maker_in_chunks_snakefile'],
            maker_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation",
            queue=config['queue'],
            jobs=config['jobs']//len(config['samples_info'])
        shell:
            """
            cd {params.maker_dir}
            snakemake -s {params.run_maker_in_chunks_snakefile} --configfile {input} --cluster "qsub -N {wildcards.sample}_chunk_annotation -q {params.queue}" -j {params.jobs} --latency-wait 60 --restart-times 3
            """


elif config['maker_parallel'] == "mpi":

    rule prep_maker_annotation_configs:
        input:
            genome=get_sample,
            liftover_gff=config["out_dir"] + "/per_sample/{sample}/MAKER_liftover/{sample}_genome_liftover.maker.output/{sample}_genome_liftover.all.gff"
        params:
            templates=config["config_templates"] + "/annotation",
            transcripts=config['annotation_transcripts_fasta'],
            proteins=config['annotation_proteins_fasta'],
            external_gff=config['external_gff']
        output:
            bopts=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_bopts.ctl",
            opts=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_opts.ctl",
            exe=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_exe.ctl"
        log:
            stdout=config["out_dir"] + "/logs/{sample}_prep_maker_annotation_configs_" + str(time()) + ".out",
            stderr=config["out_dir"] + "/logs/{sample}_prep_maker_annotation_configs_" + str(time()) + ".err"
        run:
            shell("cp {params.templates}/maker_bopts.ctl {output.bopts} > {log.stdout} 2> {log.stderr}")
            shell("cp {params.templates}/maker_exe.ctl {output.exe} >> {log.stdout} 2>> {log.stderr}")
            if params.external_gff:
                shell("sed -e 's|<genome_fasta>|{input.genome}|' -e 's|<transcripts_fasta>|{params.transcripts}|' -e 's|<proteins_fasta>|{params.proteins}|' -e 's|<gene_models_gff>|{input.liftover_gff},{params.external_gff}|' {params.templates}/maker_opts.ctl > {output.opts} 2>> {log.stderr}")
            else:
                shell("sed -e 's|<genome_fasta>|{input.genome}|' -e 's|<transcripts_fasta>|{params.transcripts}|' -e 's|<proteins_fasta>|{params.proteins}|' -e 's|<gene_models_gff>|{input.liftover_gff}|' {params.templates}/maker_opts.ctl > {output.opts} 2>> {log.stderr}")
    
    rule maker_annotation:
        input:
            bopts=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_bopts.ctl",
            opts=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_opts.ctl",
            exe=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker_exe.ctl"
        output:
            config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/annotation.done"
        log:
            index=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}_master_datastore_index.log",
            stdout=config["out_dir"] + "/logs/{sample}_maker_annotation_" + str(time()) + ".out",
            stderr=config["out_dir"] + "/logs/{sample}_maker_annotation_" + str(time()) + ".err"
        params:
            run_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation",
            base_name="{sample}",
            nodes=config['nodes'],
            ppn=config['ppn']
        threads:
            config['ppn'] * config['nodes']
    #    conda:
    #        "conda_env/maker-env.yaml"
        shell:
            """
            cd {params.run_dir}
            module load miniconda/miniconda2-4.5.4-MakerMPI
            mpirun -n {threads} -env I_MPI_FABRICS tcp maker -base {params.base_name} > {log.stdout} 2> {log.stderr}
            if [ -f {log.index} ] && [ `grep STARTED {log.index} | wc -l` == `grep FINISHED {log.index} | wc -l` ]; then touch {output}; fi
            """
    rule create_annotation_gff:
        input:
            config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}_master_datastore_index.log"
        output:
            full_gff=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.all.convert.gff",
            genes_gff=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.genes.convert.gff"
        log:
            stdout=config["out_dir"] + "/logs/{sample}_create_annotation_gff_" + str(time()) + ".out",
            stderr=config["out_dir"] + "/logs/{sample}_create_annotation_gff_" + str(time()) + ".err"
        threads:
            1
        params:
            nodes=1,
            ppn=1,
            annotation_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation"
        conda:
            "conda_env/maker-env.yaml"
        shell:
            """
            cd {params.annotation_dir}/{sample}.maker.output
            gff3_merge -d {input} -n -s > {output.full_gff} 2> {log.stderr}
            gff3_merge -d {input} -n -g -s > {output.genes_gff} 2>> {log.stderr}
            """
    
    rule create_annotation_fasta:
        input:
            config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/{sample}.maker.output/{sample}_master_datastore_index.log"
        output:
            config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.proteins.fasta"
        log:
            stdout=config["out_dir"] + "/logs/{sample}_create_annotation_fasta_" + str(time()) + ".out",
            stderr=config["out_dir"] + "/logs/{sample}_create_annotation_fasta_" + str(time()) + ".err"
        threads:
            1
        params:
            nodes=1,
            ppn=1,
            annotation_dir=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation"
        conda:
            "conda_env/maker-env.yaml"
        shell:
            """
            cd {params.annotation_dir}/{sample}.maker.output
            fasta_merge -d {input} > {log.stdout} 2> {log.stderr}
            ln {sample}.all.maker.proteins.fasta {output}
            """

rule clean_proteins_fasta:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.proteins.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.proteins.clean.fasta"
    log:
        stdout=config["out_dir"] + "/logs/{sample}_clean_proteins_fasta_" + str(time()) + ".out",
        stderr=config["out_dir"] + "/logs/{sample}_clean_proteins_fasta_" + str(time()) + ".err"
    params:
        nodes=1,
        ppn=1,
        clean_proteins_fasta_script=config['clean_proteins_fasta_script']
    conda:
        "conda_env/biopython.yaml"
    shell:
        "python {params.clean_proteins_fasta_script} {input} {output} > {log.stdout} 2> {log.stderr}"


rule proteins_busco:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.proteins.clean.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/run_BUSCO/full_table_BUSCO.tsv"
    log:
        stdout=config["out_dir"] + "/logs/{sample}_proteins_busco_" + str(time()) + ".out",
        stderr=config["out_dir"] + "/logs/{sample}_proteins_busco_" + str(time()) + ".err"
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
        python {params.busco_exe} --in {input} --out BUSCO --lineage_path {params.busco_set} --cpu {threads} --mode proteins -f > {log.stdout} 2> {log.stderr}
        """

rule proteins_interproscan:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.proteins.clean.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/interProScan/{sample}.interProScan.tsv"
    log:
        stdout=config["out_dir"] + "/logs/{sample}_proteins_interproscan_" + str(time()) + ".out",
        stderr=config["out_dir"] + "/logs/{sample}_proteins_interproscan_" + str(time()) + ".err"
    params:       
        nodes=1,
        ppn=1,
        ips_dir=config["out_dir"] + "/per_sample/{sample}/interProScan/"
    shell:
        """
        module load interproscan/5.32-71
        module load java/jdk1.8.25
        interproscan.sh -i {input} -b {params.ips_dir}/{wildcards.sample}.interProScan -t p -dp -pa -appl Pfam,ProDom,SuperFamily,PIRSF --goterms --iprlookup -f tsv,gff3 -T {params.ips_dir} > {log.stdout} 2> {log.stderr}
        """

rule proteins_blast:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.proteins.clean.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/Blast/{sample}.blast.tsv"
    log:
        stdout=config["out_dir"] + "/logs/{sample}_proteins_blast_" + str(time()) + ".out",
        stderr=config["out_dir"] + "/logs/{sample}_proteins_blast_" + str(time()) + ".err"
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
        blastp -query {input} -db {params.blast_qa_db} -evalue 1e-5 -outfmt \"6 qseqid sseqid pident qlen length mismatch gapopen evalue bitscore qcovs\" -num_threads {threads} -max_target_seqs 1 -out {output} > {log.stdout} 2> {log.stderr}
        """

rule create_repeats_gff:
    input:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.all.convert.gff"
    output:
        config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.repeats.convert.gff"
    shell:
        """
        awk '$2 == "repeatmasker" && $3 == "match"' {input} > {output}
        """

rule annotation_qa:
    input:
        genes_gff=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.genes.convert.gff",
        busco=config["out_dir"] + "/per_sample/{sample}/run_BUSCO/full_table_BUSCO.tsv",
        ips=config["out_dir"] + "/per_sample/{sample}/interProScan/{sample}.interProScan.tsv",
        blast=config["out_dir"] + "/per_sample/{sample}/Blast/{sample}.blast.tsv",
        repeats_gff=config["out_dir"] + "/per_sample/{sample}/MAKER_annotation/maker.repeats.convert.gff"
    output:
        config["out_dir"] + "/per_sample/{sample}/Annotation_QA/{sample}.QA_report.tsv"
    log:
        stdout=config["out_dir"] + "/logs/{sample}_annotation_qa_" + str(time()) + ".out",
        stderr=config["out_dir"] + "/logs/{sample}_annotation_qa_" + str(time()) + ".err"
    params:
        nodes=1,
        ppn=1,
        qa_script=config['qa_script']
    conda:
        "conda_env/annotation_qa.yaml"
    shell:
        """
        python {params.qa_script} {input.genes_gff} {output} --busco_result {input.busco} --blast_result {input.blast} --ips_result {input.ips} --repeats_gff {input.repeats_gff} > {log.stdout} 2> {log.stderr}
        """
