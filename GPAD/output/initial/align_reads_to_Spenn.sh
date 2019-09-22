#!/bin/bash
#PBS -S /bin/bash
#PBS -N BWA_align_Spenn
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/BWA_align_Spenn.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/BWA_align_Spenn.OU
#PBS -l nodes=1:ppn=10

source ~/.bashrc
hostname

cd /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial
#bwa mem -t 10 ../../data/Spenn.fasta /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/data/SRR1572628_1.fastq /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/data/SRR1572628_2.fastq > align_to_Spenn.sam
samtools view -b -h align_to_Spenn.sam > align_to_Spenn.bam
samtools sort align_to_Spenn.bam align_to_Spenn.sort -@ 10
