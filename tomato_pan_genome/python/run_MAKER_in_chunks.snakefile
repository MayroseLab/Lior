from snakemakeUtils import *
from time import time

def init(): 
    #load_info_file
    config['chunks_info'] = SampleInfoReader.sample_table_reader(filename=config['chunks_info_file'], 
                delimiter='\t', key_name='chunk', col_names=['path'])

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

localrules: all, prep_maker_configs
all_chunks = config['chunks_info'].keys()

rule all:
    input:
        config["out_dir"] + "/maker.all.convert.gff",
        config["out_dir"] + "/maker.genes.convert.gff",
        config["out_dir"] + "/maker.transcripts.fasta",
        config["out_dir"] + "/maker.proteins.fasta"
	
def get_chunk(wildcards):
    return config['chunks_info'][wildcards.chunk]['path']

rule prep_maker_configs:
    input:
        get_chunk
    output:
        bopts=config["out_dir"] + "/chunks/{chunk}/maker_bopts.ctl",
        opts=config["out_dir"] + "/chunks/{chunk}/maker_opts.ctl",
        exe=config["out_dir"] + "/chunks/{chunk}/maker_exe.ctl"
    params:
        templates=config["config_templates"],
        config_edit_script=config["config_edit_script"],
        config_kv_pairs=config["config_kv_pairs"]
    shell:
        """
        cp {params.templates}/maker_bopts.ctl {output.bopts} 
        cp {params.templates}/maker_exe.ctl {output.exe}
        python {params.config_edit_script} {params.templates}/maker_opts.ctl {output.opts} --edits genome={input} {params.config_kv_pairs}
        """

rule run_maker:
    input:
        bopts=config["out_dir"] + "/chunks/{chunk}/maker_bopts.ctl",
        opts=config["out_dir"] + "/chunks/{chunk}/maker_opts.ctl",
        exe=config["out_dir"] + "/chunks/{chunk}/maker_exe.ctl"
    output:
        config["out_dir"] + "/chunks/{chunk}/maker.done"
    log:
        index=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk_master_datastore_index.log"
    params:
        run_dir=config["out_dir"] + "/chunks/{chunk}",
    shell:
        """
        cd {params.run_dir}
        module load miniconda/miniconda2-4.5.4-MakerMPI
        maker -b chunk
        if [ -f {log.index} ] && [ `grep STARTED {log.index} | wc -l` == `grep FINISHED {log.index} | wc -l` ]; then touch {output}; fi
        """

rule create_full_gff:
    input:
       done=config["out_dir"] + "/chunks/{chunk}/maker.done",
       index=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk_master_datastore_index.log" 
    output:
        config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.gff"
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        gff3_merge -d {input.index} -n -s > {output}
        """

rule create_genes_gff:
    input:
       done=config["out_dir"] + "/chunks/{chunk}/maker.done",
       index=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk_master_datastore_index.log"
    output:
        config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.genes.gff"
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        gff3_merge -d {input.index} -n -g -s > {output}
        """

rule merge_full_gff:
    input:
       expand(config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.gff", chunk=all_chunks)
    output:
        config["out_dir"] + "/maker.all.gff"
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        gff3_merge -n -s {input} | grep -v '###' > {output}
        """

rule merge_genes_gff:
    input:
       expand(config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.genes.gff", chunk=all_chunks)
    output:
        config["out_dir"] + "/maker.genes.gff"
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        gff3_merge -n -s {input} | grep -v '###' > {output}
        """

rule convert_gff_coords:
    input:
        config["out_dir"] + "/maker.{type}.gff"
    output:
        config["out_dir"] + "/maker.{type}.convert.gff"
    params:
        coord_conversion_script = config["coord_conversion_script"]
    shell:
        "python {params.coord_conversion_script} {input} {output}"

rule create_fasta:
    input:
       done=config["out_dir"] + "/chunks/{chunk}/maker.done",
       index=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk_master_datastore_index.log",
       genes_gff=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.genes.gff"
    output:
        trans=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.maker.transcripts.fasta",
        prot=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.maker.proteins.fasta"
    params:
        out_dir = config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/",
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        cd {params.out_dir}
        fasta_merge -d {input.index}
        if [ ! -f {output.prot} ] && [ `grep -v '#' {input.genes_gff} | wc -l` == 0 ]; then touch {output.prot}; fi
        if [ ! -f {output.trans} ] && [ `grep -v '#' {input.genes_gff} | wc -l` == 0 ]; then touch {output.trans}; fi
        """

rule merge_transcripts_fasta:
    input:
       expand(config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.maker.transcripts.fasta", chunk=all_chunks)
    output:
        config["out_dir"] + "/maker.transcripts.fasta"
    shell:
        """
        cat {input} > {output}
        """
        
rule merge_proteins_fasta:
    input:
       expand(config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.maker.proteins.fasta", chunk=all_chunks)
    output:
        config["out_dir"] + "/maker.proteins.fasta"
    shell:
        """
        cat {input} > {output}
        """
