#!/bin/bash
#PBS -S /bin/bash
#PBS -N bedtools_intersect
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/bedtools_intersect.ER
#PBS -o /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial/bedtools_intersect.OU

source ~/.bashrc
hostname

cd /groups/itay_mayrose/nosnap/liorglic/Projects/GPAD/output/initial
bedtools intersect -a ITAG4.0_gene_models_exons.bed -b align_to_ref.bedGraph -wao > exons_coverage.tsv
