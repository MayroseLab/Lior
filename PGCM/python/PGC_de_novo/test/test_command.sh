#!/bin/bash
#PBS -S /bin/bash
#PBS -N PGC_de_novo_test
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/PGC_de_novo/test/PGC_de_novo_test.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/PGC_de_novo/test/PGC_de_novo_test.OU

source ~/.bashrc
hostname
conda activate snakemake

cd /groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/PGC_de_novo/test
snakefile="/groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/PGC_de_novo/PGC_de_novo.snakefile"
qsub_script="/groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/pbs_qsub_snakemake_wrapper.py"
job_script="/groups/itay_mayrose/nosnap/liorglic/Projects/PGCM/python/jobscript.sh"
snakemake -s $snakefile --configfile test.conf.yml --cluster "python $qsub_script" --latency-wait 60 --use-conda -p -j 10 --jobscript "$job_script" >out 2>err
