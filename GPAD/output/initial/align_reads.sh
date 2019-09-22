#!/bin/bash
#PBS -S /bin/bash
#PBS -N BWA_align
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/BWA_align.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/BWA_align.OU
#PBS -l nodes=1:ppn=10

source ~/.bashrc
hostname

cd /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial
bwa mem -t 10 /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/data/S_lycopersicum_chromosomes.4.00.fa /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/data/SRR1572628_1.fastq /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/data/SRR1572628_2.fastq > align_to_ref.sam
