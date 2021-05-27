#--------------------------------------------  CRP EWAS  related analysis
#---------------------------------------------------------------------------
##-------------------- overview
1. please restict to 1511 DMR loci (you can use the same lists you use to query your methQTL database)
2. residualize data (this is because we are using data after CPACOR QC and want to rid of the very long regression models arising from this)
3. We will use the residualized data as outcome for regression analysis in MR
## ------------  prerequisites
your phenotype files needs to contain 
  "log(CRP)" values
  control probes from CPACOR pipeline
  BMI
  Smoking (current, former, never)
  Houseman estimates
your DNAmethylation file needs to contain
  DNA methylation beta values from CPACOR pipeline
#-------------------------------- example residulaize script:
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
                     
write.table(res,file="nfbc66_cpacor_resids.txt",sep="\t",col.names=T,row.names=T,quote=F)
