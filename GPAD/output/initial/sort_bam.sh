#!/bin/bash
#PBS -S /bin/bash
#PBS -N sort_bam
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/sort_bam.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/sort_bam.OU
#PBS -l nodes=1:ppn=10

source ~/.bashrc
hostname

cd /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial
#samtools view -b -h align_to_ref.sam > align_to_ref.bam
samtools sort align_to_ref.bam align_to_ref.sort -@ 10
