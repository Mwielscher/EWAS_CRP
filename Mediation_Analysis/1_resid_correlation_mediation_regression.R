#--------------------------------------------  CRP EWAS  related analysis
#---------------------------------------------------------------------------
##-------------------- overview
# 1. please restict to 1511 DMR loci (you can use the same lists you use to query your methQTL database)
# 2. residualize data (this is because we are all using data after CPACOR QC and want to rid of the very long regression models arising from this)
# 3. report correlation between markers
# 4. run Regression needed for Mediaion anlysis
## ------------  prerequisites
#your phenotype files needs to contain 
#  "log(CRP)" values
#  control probes from CPACOR pipeline
#  BMI
#  Smoking (current, former, never)
#  Houseman estimates
#   your DNAmethylation file needs to contain
#  DNA methylation beta values from CPACOR pipeline
#-------------------------------- example residualize script:
lfla=as.formula('beta[i, ] ~ phe$PC1_cp + phe$PC2_cp + 
                phe$PC3_cp + phe$PC4_cp + phe$PC5_cp + phe$PC6_cp + phe$PC7_cp + 
                phe$PC8_cp + phe$PC9_cp + phe$PC10_cp + phe$PC11_cp + phe$PC12_cp + 
                phe$PC13_cp + phe$PC14_cp + phe$PC15_cp + phe$PC16_cp + phe$PC17_cp + 
                phe$PC18_cp + phe$PC19_cp + phe$PC20_cp + phe$PC21_cp + phe$PC22_cp + 
                phe$PC23_cp + phe$PC24_cp + phe$PC25_cp + phe$PC26_cp + phe$PC27_cp + 
                phe$PC28_cp + phe$PC29_cp + phe$PC30_cp')

# regression
phe=read.table("data/phe.txt",sep="\t",header=T) 
samples=as.character(phe$meth_ID)
samples=paste("X",samples,sep="")
beta=read.table("nfbc66_cpacor_betas.txt",header=T) # this assumes a subset dataset
beta=beta[,as.character(samples)]
beta=as.matrix(beta)
nvar = nrow(beta)
res=matrix(ncol=ncol(beta), nrow=nrow(beta))
rownames(res)=rownames(beta)
colnames(res)=colnames(beta)
for(i in 1:nvar) {
  tryCatch({fit = lm(lfla,na.action=na.exclude)}, error = function(error) {return(NA)})
    if(!exists("fit")){
        res[i,colnames(as.matrix(beta))] = rep(NA, length(colnames(as.matrix(beta))))
         }else{
          res[i,rownames(as.matrix(fit$residuals))] = fit$residuals
        rm(fit)
         }		
  }
                     
#write.table(res,file="nfbc66_cpacor_resids.txt",sep="\t",col.names=T,row.names=T,quote=F)

#---------------------------------------------  report correlation
# for this we just need to load the residuals and flip them
dat86=read.table("nfbc1986_age16_cpacor_resids.txt",header=T)
dat86=t(as.matrix(dat86)) 
corrNFBC86=cor(dat86,use="pairwise.complete.obs",method="pearson")
write.csv(corrNFBC86,file="CpG_correlation_NFBC1986.csv")
#---------------------------------------------------------------------------------------
#-------------------------------------------------------------   MEDIATION analysis
# data in
phe=read.table("phe_crp_NFBC66.txt",sep="\t",header=T)
meth=as.data.frame(t(read.table("nfbc66_cpacor_resids.txt",sep="\t",header=T,row.names=1)))
CpGs=as.character(colnames(meth))
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

#------------------   model 1. BMI -> CRP -> DNAmeth
# prep output files
resC=as.data.frame(matrix(NA,nrow=length(CpGs),ncol=5))
colnames(resC)=c("Estimate","Std. Error","t value","Pr(>|t|)","nsamp" )
rownames(resC)=CpGs
resA=resB=resC2=resC
i=0
for (k in CpGs){
  print(c(k))
  i=i+1
  ###  ---------------------   models
  cmod=paste(k,"~ BMI + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")            # C path
  amod=c("lnCRP ~ BMI + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos")                  # A path
  bmod=paste(k,"~ lnCRP +BMI + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")    # B path
  c2mod=paste(k,"~ BMI +lnCRP + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")    # C' path
  
  resC[i,] = tryCatch({c(summary(lm(cmod,data=dat))$coefficients[2,] ,nobs(lm(cmod,data=dat)))},error = function(error) {return(rep(NA,5))})
  resA[i,] = tryCatch({c(summary(lm(amod,data=dat))$coefficients[2,] ,nobs(lm(amod,data=dat)))},error = function(error) {return(rep(NA,5))})
  resB[i,] = tryCatch({c(summary(lm(bmod,data=dat))$coefficients[2,] ,nobs(lm(bmod,data=dat)))},error = function(error) {return(rep(NA,5))})
  resC2[i,] = tryCatch({c(summary(lm(c2mod,data=dat))$coefficients[2,] ,nobs(lm(c2mod,data=dat)))},error = function(error) {return(rep(NA,5))})
  
}
write.csv(resC,file="model_1_pathC_nfbc1966.csv")
write.csv(resC2,file="model_1_pathC2_nfbc1966.csv")
write.csv(resA,file="model_1_pathA_nfbc1966.csv")
write.csv(resB,file="model_1_pathB_nfbc1966.csv")
#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------     SMOKING MODELS
dat=dat[!(dat$smoke %in% c(2)),]  # exculde former somkers
dat$smoke[dat$smoke == 3]=2     ## code smoking to 0 and 1
dat$smoke=dat$smoke -1  
# ----------------------------------------------------------------   model 3 smoking -> CRP -> DNAmeth
# prep output files
resC=as.data.frame(matrix(NA,nrow=length(CpGs),ncol=5))
colnames(resC)=c("Estimate","Std. Error","t value","Pr(>|t|)","nsamp" )
rownames(resC)=CpGs
resA=resB=resC2=resC
i=0
for (k in CpGs){
  print(c(k))
  i=i+1
  ###  ---------------------   models
  cmod=paste(k,"~ smoke + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")       # C path
  amod=paste("lnCRP ~ smoke + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep="")         # A path
  bmod=paste(k, "~lnCRP +smoke + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")     # B path
  c2mod=paste(k,"~smoke + lnCRP + age + sex + CD4T + NK + Bcell + Mono + Neu + Eos",sep=" ")      # C' path
  
  resC[i,] = tryCatch({c(summary(lm(cmod,data=dat))$coefficients[2,] ,nobs(lm(cmod,data=dat)))},error = function(error) {return(rep(NA,5))})
  resA[i,] = tryCatch({c(summary(lm(amod,data=dat))$coefficients[2,] ,nobs(lm(amod,data=dat)))},error = function(error) {return(rep(NA,5))})
  resB[i,] = tryCatch({c(summary(lm(bmod,data=dat))$coefficients[2,] ,nobs(lm(bmod,data=dat)))},error = function(error) {return(rep(NA,5))})
  resC2[i,] = tryCatch({c(summary(lm(c2mod,data=dat))$coefficients[2,] ,nobs(lm(c2mod,data=dat)))},error = function(error) {return(rep(NA,5))})
  
}
write.csv(resC,file="model_3_pathC_nfbc1966.csv")
write.csv(resC2,file="model_3_pathC2_nfbc1966.csv")
write.csv(resA,file="model_3_pathA_nfbc1966.csv")
write.csv(resB,file="model_3_pathB_nfbc1966.csv")
