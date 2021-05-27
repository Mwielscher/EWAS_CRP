
## This script uses individual level data for 3 Cohorts (NFBC1966, NFBC1986 and Airwave). 
## For KORA we used precomputed regresssion results -- see script 1)

rm(list=ls())
#------------------------------------------------------
#------------------------  nfbc66
phe=read.table("phe_crp_NFBC66.txt",sep="\t",header=T)
meth=as.data.frame(t(read.table("nfbc66_cpacor_resids.txt",sep="\t",header=T,row.names=1)))
meth$ID=gsub("X","",as.character(rownames(meth)))
dat=merge(phe,meth,by.x="meth_ID",by.y="ID",all.x=T)
#-----------  outlier 
dat=dat[!(is.na(dat$BMI)|is.na(dat$lnCRP)), ]

table(dat$BMI>median(dat$BMI,na.rm=T)+4*(sd(dat$BMI,na.rm=T)) 
      | dat$BMI<median(dat$BMI,na.rm=T)-4*(sd(dat$BMI,na.rm=T)))
table(dat$lnCRP>median(dat$lnCRP,na.rm=T)+4*(sd(dat$lnCRP,na.rm=T)) |
        dat$lnCRP<median(dat$lnCRP,na.rm=T)-4*(sd(dat$lnCRP,na.rm=T)))

dat= dat[!(dat$BMI>median(dat$BMI,na.rm=T)+4*(sd(dat$BMI,na.rm=T)) 
           |dat$BMI<median(dat$BMI,na.rm=T)-4*(sd(dat$BMI,na.rm=T))),]
## restirct to CpGs present in all sets -- here called consensus CpGs
cons=read.table("consensus_CpGs.txt",header=F)
int1=dat[,1:28]
int2=dat[,colnames(dat) %in% cons$V1]
dat=as.data.frame(cbind(int1,int2))
#------------------------------------------------------
#------------------------  nfbc86
phe=read.table("phe_crp_NFBC86.txt",sep="\t",header=T)
colnames(phe)[4]=c("sex")
meth=as.data.frame(t(read.table("nfbc86_cpacor_resids.txt",sep="\t",header=T,row.names=1)))
meth$ID=gsub("X","",as.character(rownames(meth)))
dat86=merge(phe,meth,by.x="meth_ID",by.y="ID",all.x=T)
#-----------  outlier 
dat86=dat86[!(is.na(dat86$BMI)|is.na(dat86$lnCRP)), ]

table(dat86$BMI>median(dat86$BMI,na.rm=T)+4*(sd(dat86$BMI,na.rm=T)) 
      | dat86$BMI<median(dat86$BMI,na.rm=T)-4*(sd(dat86$BMI,na.rm=T)))
table(dat86$lnCRP>median(dat86$lnCRP,na.rm=T)+4*(sd(dat86$lnCRP,na.rm=T)) |
        dat86$lnCRP<median(dat86$lnCRP,na.rm=T)-4*(sd(dat86$lnCRP,na.rm=T)))

dat86= dat86[!(dat86$BMI>median(dat86$BMI,na.rm=T)+4*(sd(dat86$BMI,na.rm=T)) 
               |dat86$BMI<median(dat86$BMI,na.rm=T)-4*(sd(dat86$BMI,na.rm=T))),]

cons=read.table("consensus_CpGs.txt",header=F)
int1=dat86[,1:25]
int2=dat86[,colnames(dat86) %in% cons$V1]
dat86=as.data.frame(cbind(int1,int2))
#------------------------------------------------
#---------------------    AIRWAVE
phe=read.table("AIRWAVE_sens_phenotype.txt",sep="\t",header=T)
not=names(phe)[17:464]
phe=phe[,!colnames(phe) %in% not]
colnames(phe)[13]=c("sex")
meth=as.data.frame(t(read.table("AIRWAVE_CRP_CPACOR_RESIDS.txt",header=T,row.names=1)))
meth$ID=gsub("X","",as.character(rownames(meth)))
datAIR=merge(phe,meth,by.x="SAMPLE_ID",by.y="ID",all.x=T)
#-----------  outlier 
write.table(datAIR,file="AIRWAVE_CRP_EWAS_DNAmethylation.txt",col.names=T,row.names=F,quote=F,sep=" ")
datAIR=datAIR[!(is.na(datAIR$BMI)|is.na(datAIR$lnCRP)), ]

table(datAIR$BMI>median(datAIR$BMI,na.rm=T)+4*(sd(datAIR$BMI,na.rm=T)) 
      | datAIR$BMI<median(datAIR$BMI,na.rm=T)-4*(sd(datAIR$BMI,na.rm=T)))
table(datAIR$lnCRP>median(datAIR$lnCRP,na.rm=T)+4*(sd(datAIR$lnCRP,na.rm=T)) |
        datAIR$lnCRP<median(datAIR$lnCRP,na.rm=T)-4*(sd(datAIR$lnCRP,na.rm=T)))

datAIR= datAIR[!(datAIR$BMI>median(datAIR$BMI,na.rm=T)+4*(sd(datAIR$BMI,na.rm=T)) 
                 |datAIR$BMI<median(datAIR$BMI,na.rm=T)-4*(sd(datAIR$BMI,na.rm=T))),]
## restirct to CpGs present in all sets
int1=datAIR[,1:22]
int2=datAIR[,colnames(datAIR) %in% cons$V1]
datAIR=as.data.frame(cbind(int1,int2))
## --------------------------------------------------
##  --------- get KORA data
kora_pathA=read.csv("model_1_pathA_KORA_F4.csv",header=T,row.names=1)
kora_pathB=read.csv("model_1_pathB_KORA_F4.csv",header=T,row.names=1)
kora_pathC=read.csv("model_1_pathC_KORA_F4.csv",header=T,row.names=1)
kora_pathC2=read.csv("model_1_pathC2_KORA_F4.csv",header=T,row.names=1)
#-------------------------------------------------------------  start analysis
### =============================  model 1
library(meta)
methIDs=as.character(cons$V1)
output=as.data.frame(matrix(NA,nrow=c(length(methIDs)),ncol=c(20)))
rownames(output)=methIDs
colnames(output)=c("ID","pathA_coef","pathA_SE","pathA_Pval","pathB_coef","pathB_SE","pathB_P","pathC_coef",
                   "pathC_SE","pathC_P","pathC2_coef","pathC2_SE","pathC2_P","diffEstC",
                   "indirectEFF","Zscore","P_val_med","P_val_gen","total_EFF","directEFF")


for (k in methIDs){
  print(c(k))
## --------------------------  
### c- path
 cfit4=kora_pathC[rownames(kora_pathC) %in% as.character(k),]
cmod=paste(k,"~ BMI + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")   
cmod_air=paste(k,"~ BMI + age + sex + CD4T + NK + Bcell + Mono + Gran",sep=" ")  
cfit1=summary(lm(cmod,data=dat))
cfit2=summary(lm(cmod,data=dat86))
cfit3=summary(lm(cmod_air,data=datAIR))
cmeta=metagen(c(coef(cfit1)[,"Estimate"][2],coef(cfit2)[,"Estimate"][2],coef(cfit3)[,"Estimate"][2],cfit4$Estimate),
              c(coef(cfit1)[,"Std. Error"][2],coef(cfit2)[,"Std. Error"][2],coef(cfit3)[,"Std. Error"][2],cfit4$Std..Error),
              comb.random = T)
cpath=c(cmeta$TE.fixed,cmeta$seTE.fixed,cmeta$pval.fixed)
##--------------------------------
###  a path
amod=c("lnCRP ~ BMI + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos")     
amod_air=c("lnCRP ~ BMI + age + sex + CD4T + NK + Bcell + Mono + Gran")     
afit1=summary(lm(amod,data=dat))
afit2=summary(lm(amod,data=dat86))
afit3=summary(lm(amod_air,data=datAIR))
afit4=kora_pathA[rownames(kora_pathA) %in% as.character(k),]
ameta=metagen(c(coef(afit1)[,"Estimate"][2],coef(afit2)[,"Estimate"][2],coef(afit3)[,"Estimate"][2],afit4$Estimate),
              c(coef(afit1)[,"Std. Error"][2],coef(afit2)[,"Std. Error"][2],coef(afit3)[,"Std. Error"][2],cfit4$Std..Error),
              comb.random = T)
## -------------------------------
## b path
bmod=paste(k,"~ lnCRP +BMI + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")    
bmod_air=paste(k,"~ lnCRP +BMI + age + sex + CD4T + NK + Bcell + Mono + Gran",sep=" ")    
bfit1=summary(lm(bmod,data=dat))
bfit2=summary(lm(bmod,data=dat86))
bfit3=summary(lm(bmod_air,data=datAIR))
bfit4=kora_pathB[rownames(kora_pathB) %in% as.character(k),]
bmeta=metagen(c(coef(bfit1)[,"Estimate"][2],coef(bfit2)[,"Estimate"][2],coef(bfit3)[,"Estimate"][2],bfit4$Estimate),
              c(coef(bfit1)[,"Std. Error"][2],coef(bfit2)[,"Std. Error"][2],coef(bfit3)[,"Std. Error"][2],bfit4$Std..Error),
              comb.random = T)
## --------------------------------
## c'-path
c2mod=paste(k,"~ BMI +lnCRP + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")    # C' path
c2mod_air=paste(k,"~ BMI +lnCRP + age + sex + CD4T + NK + Bcell + Mono + Gran",sep=" ")
c2fit1=summary(lm(c2mod,data=dat))
c2fit2=summary(lm(c2mod,data=dat86))
c2fit3=summary(lm(c2mod_air,data=datAIR))
c2fit4=kora_pathB[rownames(kora_pathC2) %in% as.character(k),]
c2meta=metagen(c(coef(c2fit1)[,"Estimate"][2],coef(c2fit2)[,"Estimate"][2],coef(c2fit3)[,"Estimate"][2],c2fit4$Estimate),
               c(coef(c2fit1)[,"Std. Error"][2],coef(c2fit2)[,"Std. Error"][2],coef(c2fit3)[,"Std. Error"][2],c2fit4$Std..Error),
               comb.random = T)
## -----------------------------------
## Aroian Sobel  (uses Aroian test equation to calculate Z scroes)
a=ameta$TE.fixed
b=bmeta$TE.fixed
SEa=ameta$seTE.fixed
SEb=bmeta$seTE.fixed
Zscore=(a*b)/sqrt((b^2*SEa^2)+(a^2*SEb^2))+(SEa*SEb)
P_val_med=pnorm(Zscore,lower.tail=T)*2
P_val_gen=pnorm(-abs(Zscore),lower.tail=T)*2
total=cmeta$TE.fixed
direct=c2meta$TE.fixed
indirect=a*b
###  --------------  populate output file
output[k,"ID"]=c(k)
output[k,2:7]=c(a,SEa,ameta$pval.fixed,b,SEb,bmeta$pval.fixed)
output[k,8:14]=c(cpath, 
                 c2meta$TE.fixed,c2meta$seTE.fixed,c2meta$pval.fixed,
                 cmeta$TE.fixed - c2meta$TE.fixed)
output[k,15:ncol(output)]=c(indirect,Zscore,P_val_med,P_val_gen,total,direct)

rm(Zscore,a,b,SEa,SEb,indirect,total,direct)
}

write.csv(output,file="meta_analysed_BMI_via_CRP_to_DNAmeth_mediation_FINAL_all_cohorts.csv")

## ------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#----------------------                          model 3 smoking -> crp -> DNAmeth

dat=dat[!(dat$smoke %in% c(2)),]  # exculde former somkers
dat$smoke[dat$smoke == 3]=2
dat$smoke=dat$smoke -1
####
dat86=dat86[!(dat86$smoke %in% c(2)),]  # exculde former somkers
dat86$smoke[dat86$smoke == 3]=2
dat86$smoke=dat86$smoke -1
###
colnames(datAIR)[22]=c("smoke")
datAIR=datAIR[!(datAIR$smoke %in% c(2)),]  # exculde former somkers
datAIR$smoke[datAIR$smoke == 3]=2
datAIR$smoke=datAIR$smoke -1
###############
##  --------- get KORA
kora_pathA=read.csv("model_3_pathA_KORA_F4.csv",header=T,row.names=1)
kora_pathB=read.csv("model_3_pathB_KORA_F4.csv",header=T,row.names=1)
kora_pathC=read.csv("model_3_pathC_KORA_F4.csv",header=T,row.names=1)
kora_pathC2=read.csv("model_3_pathC2_KORA_F4.csv",header=T,row.names=1)

#-------------------------------------------------------------  start analysis
library(meta)
methIDs=as.character(cons$V1)
output=as.data.frame(matrix(NA,nrow=c(length(methIDs)),ncol=c(20)))
rownames(output)=methIDs
colnames(output)=c("ID","pathA_coef","pathA_SE","pathA_Pval","pathB_coef","pathB_SE","pathB_P","pathC_coef",
                   "pathC_SE","pathC_P","pathC2_coef","pathC2_SE","pathC2_P","diffEstC",
                   "indirectEFF","Zscore","P_val_med","P_val_gen","total_EFF","directEFF")


for (k in methIDs){
  print(c(k))
  ### c- path
  cmod=paste(k,"~ smoke + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")   
  cmod_air=paste(k,"~ smoke + age + sex + CD4T + NK + Bcell + Mono + Gran",sep=" ")  
  cfit1=summary(lm(cmod,data=dat))
  cfit2=summary(lm(cmod,data=dat86))
  cfit3=summary(lm(cmod_air,data=datAIR))
  cfit4=kora_pathC[rownames(kora_pathC) %in% as.character(k),]
  cmeta=metagen(c(coef(cfit1)[,"Estimate"][2],coef(cfit2)[,"Estimate"][2],coef(cfit3)[,"Estimate"][2],cfit4$Estimate),
                c(coef(cfit1)[,"Std. Error"][2],coef(cfit2)[,"Std. Error"][2],coef(cfit3)[,"Std. Error"][2],cfit4$Std..Error),
                comb.random = T)
  cpath=c(cmeta$TE.fixed,cmeta$seTE.fixed,cmeta$pval.fixed)
  ###  a path
  amod=paste("lnCRP ~ smoke + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep="")     
  amod_air=paste("lnCRP ~ smoke + age + sex + CD4T + NK + Bcell + Mono + Gran",sep="")     
  afit1=summary(lm(amod,data=dat))
  afit2=summary(lm(amod,data=dat86))
  afit3=summary(lm(amod_air,data=datAIR))
  afit4=kora_pathA[rownames(kora_pathA) %in% as.character(k),]
  ameta=metagen(c(coef(afit1)[,"Estimate"][2],coef(afit2)[,"Estimate"][2],coef(afit3)[,"Estimate"][2],afit4$Estimate),
                c(coef(afit1)[,"Std. Error"][2],coef(afit2)[,"Std. Error"][2],coef(afit3)[,"Std. Error"][2],cfit4$Std..Error),
                comb.random = T)
  
  ## b path
  bmod=paste(k, "~lnCRP +smoke + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")    
  bmod_air=paste(k, "~lnCRP +smoke + age + sex + CD4T + NK + Bcell + Mono + Gran",sep=" ")    
  bfit1=summary(lm(bmod,data=dat))
  bfit2=summary(lm(bmod,data=dat86))
  bfit3=summary(lm(bmod_air,data=datAIR))
  bfit4=kora_pathB[rownames(kora_pathB) %in% as.character(k),]
  bmeta=metagen(c(coef(bfit1)[,"Estimate"][2],coef(bfit2)[,"Estimate"][2],coef(bfit3)[,"Estimate"][2],bfit4$Estimate),
                c(coef(bfit1)[,"Std. Error"][2],coef(bfit2)[,"Std. Error"][2],coef(bfit3)[,"Std. Error"][2],bfit4$Std..Error),
                comb.random = T)
  
  ## c'-path
  c2mod=paste(k,"~smoke + lnCRP + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")    # C' path
  c2mod_air=paste(k, " ~ smoke + lnCRP + age + sex + CD4T + NK + Bcell + Mono + Gran",sep=" ")
  c2fit1=summary(lm(c2mod,data=dat))
  c2fit2=summary(lm(c2mod,data=dat86))
  c2fit3=summary(lm(c2mod_air,data=datAIR))
  c2fit4=kora_pathC2[rownames(kora_pathC2) %in% as.character(k),]
  c2meta=metagen(c(coef(c2fit1)[,"Estimate"][2],coef(c2fit2)[,"Estimate"][2],coef(c2fit3)[,"Estimate"][2],c2fit4$Estimate),
                 c(coef(c2fit1)[,"Std. Error"][2],coef(c2fit2)[,"Std. Error"][2],coef(c2fit3)[,"Std. Error"][2],c2fit4$Std..Error),
                 comb.random = T)
  
  ## Aroian Sobel  (uses Aroian test equation to calculate Z scroes)
  a=ameta$TE.fixed
  b=bmeta$TE.fixed
  SEa=ameta$seTE.fixed
  SEb=bmeta$seTE.fixed
  Zscore=(a*b)/sqrt((b^2*SEa^2)+(a^2*SEb^2))+(SEa*SEb)
  P_val_med=pnorm(Zscore,lower.tail=T)*2
  P_val_gen=pnorm(-abs(Zscore),lower.tail=T)*2
  total=cmeta$TE.fixed
  direct=c2meta$TE.fixed
  indirect=a*b
  ###  --------------  populate output file
  output[k,"ID"]=c(k)
  output[k,2:7]=c(a,SEa,ameta$pval.fixed,b,SEb,bmeta$pval.fixed)
  output[k,8:14]=c(cpath, 
                   c2meta$TE.fixed,c2meta$seTE.fixed,c2meta$pval.fixed,
                   cmeta$TE.fixed - c2meta$TE.fixed)
  output[k,15:ncol(output)]=c(indirect,Zscore,P_val_med,P_val_gen,total,direct)
  
  rm(Zscore,a,b,SEa,SEb,indirect,total,direct)
}


write.csv(output,file="meta_analysed_smoking_via_CRP_to_DNA_meth_mediation_FINAL_all_COHORTS.csv")

