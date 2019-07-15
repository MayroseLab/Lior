#!/bin/bash
#PBS -S /bin/bash
#PBS -N SNPhylo_iris
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose/liorglic/SNPhylo/SNPhylo_iris.ER
#PBS -o /groups/itay_mayrose/liorglic/SNPhylo/SNPhylo_iris.OU

source ~/.bashrc
hostname
conda activate SNPhylo
cd /groups/itay_mayrose/liorglic/SNPhylo
/groups/itay_mayrose/liorglic/SNPhylo/SNPhylo/snphylo.sh -v RLZ18_aligned_genotypes_LowMinPop.vcf -A -b -a 123427 -p 90 -M 0.5 -m 0.01 -o Iris_mesop >out 2>err
