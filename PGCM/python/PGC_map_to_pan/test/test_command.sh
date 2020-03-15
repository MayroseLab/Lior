#!/bin/bash
#PBS -S /bin/bash
#PBS -N PGC_map_to_pan_test
#PBS -r y
#PBS -q itay_25_1
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/PGC_map_to_pan/test/PGC_map_to_pan_test.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/PGC_map_to_pan/test/PGC_map_to_pan_test.OU

source ~/.bashrc
hostname
conda activate snakemake
PATH="/groups/itay_mayrose/liorglic/miniconda3/envs/snakemake/bin:$PATH"

cd /groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/PGC_map_to_pan/test
snakefile="/groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/PGC_map_to_pan/PGC_map_to_pan.snakefile"
qsub_script="/groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/pbs_qsub_snakemake_wrapper.py"
job_script="/groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/jobscript.sh"
snakemake -s $snakefile --configfile test.conf.yml --cluster "python $qsub_script" --latency-wait 60 --use-conda -p -j 10 --jobscript "$job_script" >out 2>err
