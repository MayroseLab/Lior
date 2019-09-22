#!/bin/bash
#PBS -S /bin/bash
#PBS -N BWA_index
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/BWA_index.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/BWA_index.OU

source ~/.bashrc
hostname

bwa index ../../data/Spenn.fasta
bwa index ../../data/Solanum_chilense.scaffolds.fa
bwa index ../../data/Spimp.fasta
