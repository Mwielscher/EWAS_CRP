#!/bin/bash
#PBS -lselect=1:ncpus=1:mem=5GB
#PBS -lwalltime=2:00:00
#PBS -N run_direct_66

#subimt qsub -J 1-1000 scriptname.sh
module load rvtests/2016-05-04

tid=$PBS_ARRAY_INDEX

indir='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/nfbc66/data/'
outdir='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/nfbc66/CpG_cause/raw/'

i=$(sed -n ${tid}p /rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/CpGcauseCRP_GWAS_list.txt | awk '{print $1}')
k=$(sed -n ${tid}p /rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/CpGcauseCRP_GWAS_list.txt | awk '{print $2}')

rvtest --inVcf ${indir}nfbc1966_for_meth.vcf.gz --pheno ${indir}MR_pheno_nfbc1966.pheno --pheno-name lnCRP\ 
 --covar ${indir}MR_pheno_nfbc1966.covar --covar-name ${i},sex,age,CD8T,CD4T,NK,Bcell,Mono,Neu,Eos --rangeList ${k}\
 --out ${outdir}DIRECT_EFF_${i} --kinship ${indir}nfbc1966_kinship.kinship --meta dominant --freqLower 0.025

gzip -d ${outdir}DIRECT_EFF_${i}.MetaDominant.assoc.gz
sed '/^#/ d' ${outdir}DIRECT_EFF_${i}.MetaDominant.assoc > ${outdir}DIRECT_EFF_${i}.MetaDominant
header=$(sed -n 1p ${outdir}DIRECT_EFF_${i}.MetaDominant)
sed -i '/NA/ d' ${outdir}DIRECT_EFF_${i}.MetaDominant
sed -i '1d' ${outdir}DIRECT_EFF_${i}.MetaDominant
awk '{print $0"\t"1/$14}' ${outdir}DIRECT_EFF_${i}.MetaDominant > ${outdir}DIRECT_EFF_${i}.final
echo -e "$header\tSE" | cat - ${outdir}DIRECT_EFF_${i}.final > ${outdir}DIRECT_EFF_${i}_nfbc1966_1 
awk '{print $0"\t"$1":"$2}' ${outdir}DIRECT_EFF_${i}_nfbc1966_1 > ${outdir}DIRECT_EFF_${i}_nfbc1966
rm ${outdir}DIRECT_EFF_${i}_nfbc1966_1
