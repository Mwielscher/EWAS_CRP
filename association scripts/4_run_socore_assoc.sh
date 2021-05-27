#!/bin/bash

module load R
data='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/nfbc66/data/'
wdr='/rds/general/user/mwielsch/home/WORK/CRP/mendelian_randomization/scores/CRP_cause/nfbc66/'


R --vanilla << EOF
dat=read.table("${data}MR_pheno_nfbc1966.pheno",header=T)

grs=read.table("${wdr}crp_score_nfbc66.profile",header=T)
grs1=as.data.frame(cbind(as.character(grs\$FID),grs\$SCORE))
colnames(grs1)=c("ID","CRP_score")
dat1=merge(dat,grs1,by.x="fid",by.y="ID",all.x=T)
dat1\$CRP_score=as.numeric(as.character(dat1\$CRP_score))
crp=read.table("${wdr}nfbc66_lnCRP.txt",header=T)
dat1=merge(dat1,crp,by.x="fid",by.y="ID",all.x=T)
dat1\$lnCRP=as.numeric(as.character(dat1\$lnCRP))
CpGs=read.table("${wdr}CpG_list")
outc=c(as.character(CpGs\$V1),"lnCRP")

####-----------------------   regressions
res=as.data.frame(matrix(NA,nrow=length(outc),ncol=5))
colnames(res)=c("Estimate","Std. Error","t value","P_value","nsamp" )
rownames(res)=as.character(outc)
i=0
for (k in outc){
  print(c(k))
  i=i+1
  ###  ---------------------   model
  cmod=paste(k,"~ CRP_score + sex + CD4T + NK + Bcell + Mono + Neu + Eos + PC1 + PC2 + PC3 + PC4 + PC5",sep=" ")            
  res[i,] = tryCatch({c(summary(lm(cmod,data=dat1))\$coefficients[2,] ,nobs(lm(cmod,data=dat1)))},error = function(error) {return(rep(NA,5))})
}
res=res[!is.na(res\$Estimate),]
write.table(res,file="${wdr}CRP_causes_DNAmeth_nfbc1966.txt",sep="\t",col.names=T,row.names=T,quote=F)

EOF
