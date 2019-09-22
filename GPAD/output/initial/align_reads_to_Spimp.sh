#!/bin/bash
#PBS -S /bin/bash
#PBS -N BWA_align_Spimp
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/BWA_align_Spimp.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/BWA_align_Spimp.OU
#PBS -l nodes=1:ppn=10

source ~/.bashrc
hostname

cd /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial
#bwa mem -t 10 ../../data/Spimp.fasta /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/data/SRR1572628_1.fastq /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/data/SRR1572628_2.fastq > align_to_Spimp.sam
samtools view -b -h align_to_Spimp.sam > align_to_Spimp.bam
samtools sort align_to_Spimp.bam align_to_Spimp.sort -@ 10

