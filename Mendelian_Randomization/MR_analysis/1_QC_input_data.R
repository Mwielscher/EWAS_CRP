## compare effect sizes between cohorts:
## the files are the result files from script 4_run_score_assoc.sh
rm(list=ls())
kora=read.table("CRP_causes_DNAmeth_KORA.txt",sep="\t",header=T)
kora=kora[,c("ID","Estimate")]
colnames(kora)[2]=c('kora')
nfbc66=read.table("CRP_causes_DNAmeth_nfbc1966.txt",sep="\t",header=T)
nfbc66=nfbc66[,c("ID","Estimate")]
colnames(nfbc66)[2]=c('nfbc66')
dat=merge(kora,nfbc66,by.x="ID",by.y="ID")
dat=dat[!dat$ID %in% c("lnCRP"),]
# --------------------
nfbc86=read.table("CRP_causes_DNAmeth_nfbc1986.txt",sep="\t",header=T)
nfbc86=nfbc86[,c("ID","Estimate")]
colnames(nfbc86)[2]=c('nfbc86')
dat=merge(dat,nfbc86,by.x="ID",by.y="ID")
#
airwave=read.table("CRP_causes_DNAmeth_airwave.txt",sep="\t",header=T)
airwave=airwave[,c("ID","Estimate")]
colnames(airwave)[2]=c('airwave')
dat=merge(dat,airwave,by.x="ID",by.y="ID")
#
RS=read.table("CRP_causes_DNAmeth_RS.txt",sep="\t",header=T)
RS=RS[,c("ID","Estimate")]
colnames(RS)[2]=c('RS')
dat=merge(dat,RS,by.x="ID",by.y="ID")
#
PAN=read.table("CRP_causes_DNAmeth_PAN.txt",sep="\t",header=T)
PAN=PAN[,c("ID","Estimate")]
colnames(PAN)[2]=c('PAN')
dat=merge(dat,PAN,by.x="ID",by.y="ID")
#
NTR=read.table("CRP_causes_DNAmeth_NTR.txt",sep="\t",header=T)
NTR=NTR[,c("ID","Estimate")]
colnames(NTR)[2]=c('NTR')
dat=merge(dat,NTR,by.x="ID",by.y="ID")
##
LLS=read.table("CRP_causes_DNAmeth_LLS.txt",sep="\t",header=T)
LLS=LLS[,c("ID","Estimate")]
colnames(LLS)[2]=c('LLS')
dat=merge(dat,LLS,by.x="ID",by.y="ID")
##
LL=read.table("CRP_causes_DNAmeth_LL.txt",sep="\t",header=T)
LL=LL[,c("ID","Estimate")]
colnames(LL)[2]=c('LL')
dat=merge(dat,LL,by.x="ID",by.y="ID")
##
CODAM=read.table("CRP_causes_DNAmeth_CODAM.txt",sep="\t",header=T)
CODAM=CODAM[,c("ID","Estimate")]
colnames(CODAM)[2]=c('CODAM')
dat=merge(dat,CODAM,by.x="ID",by.y="ID")
##------------  done
ref=read.csv("MR_BMI_EWAS_nature.csv",header=T)
table(dat$ID %in% ref$CpG)

library(corrplot)
#test=cor(dat1[,2:ncol(dat1)])
#corrplot(test)
test1=cor(dat[,2:ncol(dat)])
#test1[test1>0.2]=0.2
corrplot(test1,method="number",is.corr = FALSE, type="upper",tl.cex=1.5,diag=F)

library(ggsci)
library(reshape2)
library(ggplot2)
dat23<- melt(dat[,2:ncol(dat)])
colnames(dat23)=c("COHORT","coefficients")
ggplot(dat23,aes(x=coefficients, fill=COHORT)) + geom_density(alpha=0.25) + xlim(-3,3) +
  theme_bw() +
  theme(axis.text.x = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),
        axis.title.x =element_text(size=17),
        axis.title.y =element_text(size=17)) 

###  -----------------  looks ok
















