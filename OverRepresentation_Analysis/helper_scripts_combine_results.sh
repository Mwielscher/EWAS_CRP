#!/bin/bash

dir='/rds/general/user/mwielsch/home/WORK/CRP/enrichment/ALL_loci/result/'
out='/rds/general/user/mwielsch/home/WORK/CRP/enrichment/ALL_loci/'
traits='CGIcorr_POS CGIcorr_DNAse H3K27 H3K4 HiC CHROM_state encodeTF GWAS_cat_'

 
for k in ${traits}; do

echo "$k"

echo "trait	expected	expected_PERC	observed	observed_PERC	ePval_incr	ePval_decr	fPval" > ${out}FINAL_ALL_loci_res_${k}.txt
awk 'FNR == 2' ${dir}res_${k}*.txt >> ${out}FINAL_ALL_loci_res_${k}.txt


done


