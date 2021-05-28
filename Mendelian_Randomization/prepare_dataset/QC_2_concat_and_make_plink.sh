#!/bin/bash
#PBS -lselect=1:ncpus=1:mem=20GB
#PBS -lwalltime=12:00:00
#PBS -N tesing_66

module load bcftools/2015-02-17
module load htslib/1.3.2
module load plink
out='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/nfbc66/data/'


bcftools concat -f ${out}concat_input -Oz -o ${out}nfbc1966_for_meth.vcf.gz
bgzip -d ${out}nfbc1966_for_meth.vcf.gz
bgzip ${out}nfbc1966_for_meth.vcf
tabix ${out}nfbc1966_for_meth.vcf.gz

## the sort is not necessary 
#bcftools sort ${out}big_vcf.bcf -T ${out}sorting/ --max-mem 20G -Ov -o ${out}nfbc1966_for_meth_inter.vcf
#(grep ^"#" nfbc1966_for_meth_inter.vcf; grep -v ^"#" nfbc1966_for_meth_inter.vcf | sed 's:^chr::ig') | bgzip -c > ${out}nfbc1966_for_meth.vcf.gz


plink --vcf ${out}nfbc1966_for_meth.vcf.gz --make-bed --double-id --allow-extra-chr -out ${out}nfbc1966_for_meth
