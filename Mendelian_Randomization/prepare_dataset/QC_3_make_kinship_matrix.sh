#!/bin/bash
#PBS -lselect=1:ncpus=1:mem=20GB
#PBS -lwalltime=22:00:00
#PBS -N kinship_66

module load bcftools/2015-02-17
module load rvtests/2016-05-04

out='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/nfbc66/data/'


vcf2kinship --inVcf ${out}nfbc1966_for_meth.vcf.gz --minMAF 0.05 --ped ${out}nfbc1966.ped --bn --out ${out}nfbc1966_kinship
