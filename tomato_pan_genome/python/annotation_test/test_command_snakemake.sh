snakemake -s ../annotate_genome.snakefile --configfile snakemake_config.yaml --cluster "qsub -N maker_try -p 0 -q itaym -l nodes={params.nodes}:ppn={params.ppn}" -j 4 --latency-wait 60 --use-conda
