#!/bin/bash
#PBS -S /bin/bash
#PBS -N BWA_align_Schil
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/BWA_align_Schil.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/BWA_align_Schil.OU
#PBS -l nodes=1:ppn=10

source ~/.bashrc
hostname

cd /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial
#bwa mem -t 10 ../../data/Solanum_chilense.scaffolds.fa /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/data/SRR1572628_1.fastq /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/data/SRR1572628_2.fastq > align_to_Schil.sam
samtools view -b -h align_to_Schil.sam > align_to_Schil.bam
samtools sort align_to_Schil.bam align_to_Schil.sort -@ 10
