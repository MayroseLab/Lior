#!/bin/bash
#PBS -S /bin/bash
#PBS -N variant_calling_pipeline_test
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/gene_content_evolution/python/variant_calling_pipeline_test/variant_calling_pipeline_test.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/gene_content_evolution/python/variant_calling_pipeline_test/variant_calling_pipeline_test.OU

source ~/.bashrc
hostname
conda activate snakemake-try

cd /groups/itay_mayrose/nosnap/liorglic/Projects/gene_content_evolution/python/variant_calling_pipeline_test
snakemake -s ../variant_calling_pipeline.snakefile --configfile test.conf.yml --cluster "python /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/python/pbs_qsub_snakemake_wrapper.py" --latency-wait 60 --use-conda -p -j 10 --jobscript "/groups/itay_mayrose/nosnap/liorglic/Projects/tomato_pan_genome/python/jobscript.sh"
