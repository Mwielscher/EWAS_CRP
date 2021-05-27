#!/bin/bash
#PBS -lselect=1:ncpus=1:mem=5GB
#PBS -lwalltime=2:00:00
#PBS -N run_assoc_66

#subit qsub -J 1-1511 4_assoc_CpG_cause_CRP.sh


module load rvtests/2016-05-04
module load htslib/1.3.2

tid=$PBS_ARRAY_INDEX

indir='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/nfbc66/data/'
interdir='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/nfbc66/CRP_cause/raw/'
outdir='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/nfbc66/CRP_cause/'
cohort='nfbc1966'
sites='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/CRP_SNP_range_list'

i=$(sed -n ${tid}p /rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/CpGcauseCRP_GWAS_list.txt | awk '{print $1}')


rvtest --inVcf ${indir}nfbc1966_for_meth.vcf.gz --pheno ${indir}MR_pheno_nfbc1966.pheno --pheno-name ${i} --covar ${indir}MR_pheno_nfbc1966.covar --covar-name lnCRP,sex,age,CD8T,CD4T,NK,Bcell,Mono,Neu,Eos --rangeFile ${sites} --out ${interdir}DIRECT_EFF_${i} --kinship ${indir}nfbc1966_kinship.kinship --meta dominant --freqLower 0.025

gzip -d ${interdir}DIRECT_EFF_${i}.MetaDominant.assoc.gz
sed '/^#/ d' ${interdir}DIRECT_EFF_${i}.MetaDominant.assoc > ${interdir}DIRECT_EFF_${i}.MetaDominant
header=$(sed -n 1p ${interdir}DIRECT_EFF_${i}.MetaDominant)
sed -i '/NA/ d' ${interdir}DIRECT_EFF_${i}.MetaDominant
sed -i '1d' ${interdir}DIRECT_EFF_${i}.MetaDominant
awk '{print $0"\t"1/$14}' ${interdir}DIRECT_EFF_${i}.MetaDominant > ${interdir}DIRECT_EFF_${i}.final
echo -e "$header\tSE" | cat - ${interdir}DIRECT_EFF_${i}.final > ${outdir}DIRECT_EFF_${i}_${cohort}_1
awk '{print $0"\t"$1":"$2}' ${outdir}DIRECT_EFF_${i}_${cohort}_1 > ${outdir}DIRECT_EFF_${i}_${cohort}
rm ${outdir}DIRECT_EFF_${i}_${cohort}_1
 
