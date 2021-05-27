#!/bin/bash


datdir='/project/eph/CAPICE_p67795/HRCimputation/66/imputed20161005/polymorph/'
out='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/nfbc66/'

for chrom in $(seq 22) 
do
echo "this is chrom${chrom} now !!"

cat << EOF > run_chr${chrom}.sh 

#!/usr/bin/bash
#PBS -lselect=1:ncpus=1:mem=10GB
#PBS -lwalltime=08:00:00
#PBS -N testing_66

module load bcftools/2015-02-17
module load vcftools/0.1.13
module load htslib/1.3.2


bcftools view -S ${out}NFBC66_meth_samples --force-samples ${datdir}nfbc66_chr${chrom}_hrc.imputed.poly.vcf.gz -Oz -o ${out}data/final_chr${chrom}.vcf.gz

#vcftools --gzvcf ${out}data/inter_chr${chrom}.vcf.gz --mac 1 --recode --recode-INFO-all --stdout | gzip -c > ${out}data/final_chr${chrom}.vcf.gz
bgzip -d ${out}data/final_chr${chrom}.vcf.gz

bgzip ${out}data/final_chr${chrom}.vcf

tabix ${out}data/final_chr${chrom}.vcf.gz
rm ${out}data/inter_chr${chrom}.vcf.gz

EOF
chmod 755 run_chr${chrom}.sh
qsub run_chr${chrom}.sh
done

