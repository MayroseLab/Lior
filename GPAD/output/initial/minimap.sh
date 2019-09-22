#!/bin/bash
#PBS -S /bin/bash
#PBS -N minimap
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/minimap.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/minimap.OU

source ~/.bashrc
hostname
conda activate minimap2

cd /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial
minimap2 -ax splice:hq -uf ../../data/S_lycopersicum_chromosomes.4.00.fa relative_transcripts.fasta > relative_transcripts_vs_heinz.sam
