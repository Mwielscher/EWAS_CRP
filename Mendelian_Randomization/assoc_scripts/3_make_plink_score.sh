#!/bin/bash

module load plink
module load R
data='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/nfbc66/data/'
wdr='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/scores/CRP_cause/'

R --vanilla << EOF
library(data.table)
dat=fread("${data}nfbc1966_for_meth.bim")
grs=read.table("${wdr}NO_DIRECT_SNPs_initial_CRP_GWAS_score.txt")
dat1=merge(grs,dat,by.x="V1",by.y="V2",all.x=T)
colnames(dat1)=c("ID","grs_A1","grs_A2","grs_beta","grs_P","bim_chr","bim_morgans","bim_pos","bim_A1","bim_A2")
dat1\$grs_A1=as.character(dat1\$grs_A1)
dat1\$grs_A2=as.character(dat1\$grs_A2)
dat1\$bim_A1=as.character(dat1\$bim_A1)
dat1\$bim_A2=as.character(dat1\$bim_A2)

for (k in 1:nrow(dat1)){
	if (toupper(dat1\$grs_A1[k]) %in% dat1\$bim_A1[k]) {
	dat1\$grs_A1[k]=dat1\$bim_A1[k]
	dat1\$grs_A2[k]=dat1\$bim_A2[k]
	}else {
	dat1\$grs_A1[k]=dat1\$bim_A1[k]
	dat1\$grs_A2[k]=dat1\$bim_A2[k]
	dat1\$grs_beta[k]=as.numeric(as.character(dat1\$grs_beta[k]))*-1
	}
}
###  you can also let plink do the allele alingment for you and skipt the R part of the script

ID_new=paste(dat1\$ID,"_",dat1\$grs_A1,"_",dat1\$grs_A2,sep="")

out=as.data.frame(cbind(ID_new,dat1\$grs_A1,dat1\$grs_beta))

write.table(out,file="${wdr}CRP_GWAS_for_nfbc1966.grs",sep="\t",col.names=F,row.names=F,quote=F)
EOF

awk '{print $1"\t"$2"_"$5"_"$6"\t"$3"\t"$4"\t"$5"\t"$6}' ${data}nfbc1966_for_meth.bim > ${data}nfbc1966_for_meth_2.bim

plink --bed ${data}nfbc1966_for_meth.bed --bim ${data}nfbc1966_for_meth_2.bim --fam ${data}nfbc1966_for_meth.fam --score ${wdr}CRP_GWAS_for_nfbc1966.grs --out ${wdr}crp_score_nfbc66 
