#!/bin/bash
#PBS -S /bin/bash
#PBS -N make_bedGraph 
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/make_bedGraph.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/make_bedGraph.OU

source ~/.bashrc
hostname

cd /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial
bedtools genomecov -bga -split -ibam align_to_Spimp.sort.bam > align_to_Spimp.bedGraph
bedtools genomecov -bga -split -ibam align_to_Spenn.sort.bam > align_to_Spenn.bedGraph
bedtools genomecov -bga -split -ibam align_to_Schil.sort.bam > align_to_Schil.bedGraph
