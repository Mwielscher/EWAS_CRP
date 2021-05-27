rm(list=ls())
####   read in list of CRP assoiated CpGs from EWAS
tt.crp=read.csv("marker_master_list_may2020.csv")
#------------------------------------------------------------------------------------------
#------------------------------------------                   data in   
load("NFBC66_resid_age46_CRP_hits.RData")  # batch - sex - age - houseman corrected methylation data from NFBC1966 age 46:
##  ------- I am using residulised DNA methylation from an EPIC array (DNAmethylation beta values, produced with CPACOR pipeline)
##  --------- I attached a few lines of code that do the residualisation of the data to the end of this script - 
## ------------- This way you can remove the infuence sex, age and cell type distribution might have on the data
#-- If that does work with your data or is already done in some other way, you can also use normalised DNA methylaiton beta values, because the data will be scaled anyway
# -------  The idea is that we do not have to add covariates to the risk score association model.
phe=read.csv("phe_T2.csv")                # phenotype file NFBC1966
# ----------  there are example files in the zipped folder I send
# --------------------------------------            align input data
dat.t2=resid_46[rownames(resid_46) %in% as.character(tt.crp$ID),]   #.t2 because those are data from the 2nd data colletion
phe=phe[phe$meth_ID %in% colnames(dat.t2),]
tt.crp=tt.crp[as.character(tt.crp$ID) %in% row.names(dat.t2),]
rownames(dat.t2)=as.character(tt.crp$ID)           
dat.t2=dat.t2[,as.character(phe$meth_ID)]
# -----------------------------------------------------------------------------------------
#-------------------------------------------       loop through marker subsets and clinical outcomes
outcomes=c("COPD","malign_neoplasm","T2D")
# ----------------------------------------------------------
subsets=c("all_loci","CpG_mediated_BMI","CpG_mediated_smoke","CORRELATION_cluster_1","CORRELATION_cluster_2","CpG_CAUSE_CRP","GeneExpr")
res=as.data.frame(matrix(NA,nrow=(length(subsets)*length(outcomes)),ncol=10))
colnames(res)=c("ID", "trait","markerList", "Estimate","Std. Error","t value","Pr(>|t|)","nsamp","score_total","scoreCpGs" )
y=1
##   --------------------------------------   subset loci loop starts here
for (u in subsets){
  print(paste0("make CpG score for ",u))
  marker=as.character(tt.crp$ID[tt.crp[,u] %in% c(1)])
  tt=tt.crp[tt.crp$ID %in% marker,]
  dat.t2=resid_46[rownames(resid_46) %in% marker,]   #.t2 because those are data from the 2nd data collection
  phe=phe[phe$meth_ID %in% colnames(dat.t2),]
  tt=tt[as.character(tt$ID) %in% row.names(dat.t2),]
  rownames(dat.t2)=as.character(tt$ID)           
  dat.t2=dat.t2[,as.character(phe$meth_ID)]
  score.l=nrow(dat.t2)
  score.cpg=paste(as.character(rownames(dat.t2)),collapse=":")
  # ----------------------------------------------------------------------------------------
  # -------------------------------------  scale DNA methylation data
  scaled_dat.t2 <- matrix(0, nrow(dat.t2), ncol(dat.t2))
  for (d in 1:nrow(dat.t2)){
      scaled_dat.t2[d,] <- scale(dat.t2[d,]) 
    }
  rownames(scaled_dat.t2) <- rownames(dat.t2)
  colnames(scaled_dat.t2) <- colnames(dat.t2)
  #--------------------------------       calculate meth score
    scaled_dat.t2[is.na(scaled_dat.t2)]=0
    meth_score <- NULL
    for (i in 1:ncol(scaled_dat.t2)) {
        meth_score[i] <- as.numeric(scaled_dat.t2[,i]) %*% tt$TE_Effect  ## set of marker coeffs
    }
    scaled_meth_score <- as.data.frame(scale(meth_score))
    rownames(scaled_meth_score)=colnames(scaled_dat.t2)
    scaled_meth_score$ID=rownames(scaled_meth_score)
    colnames(scaled_meth_score)=c("SCORE","ID")
  # -----------------------------------------------------------------------------------------
  #-------------------------------------------       association testing
   for (k in outcomes){ 
      scaled_meth_score=scaled_meth_score[as.character(phe$meth_ID),]  ## -- should be aligned already
      res[y,1:3]=c(paste(k,u,sep="_"),as.character(k),as.character(u))
      res[y,4:8] = tryCatch({c(summary(glm(phe[,k] ~ scaled_meth_score$SCORE,family=binomial))$coefficients[2,] 
                                   ,nobs(glm(phe[,k] ~ scaled_meth_score$SCORE,family=binomial)))},error = function(error) {return(rep(NA,5))})
      res[y,9:10]=c(score.l,score.cpg)
      y=y+1
    } #  outcomes
    rm(score.l,score.cpg)
} # marker subset loop
# --------  make result file
write.csv(res,file="yourCohort_yourInitials_date_CRP_risk_score.csv",row.names=F)
#----------------------------------------------    END
#-----------------------------------------------------------------------------------------------------------------


#-------------------------------- example residualize script:
lfla=as.formula('beta[i, ] ~ phe$age + phe$sex + phe$CD4T + phe$NK + phe$Bcell + phe$Mono + phe$Neu + phe$Eos +
                phe$PC1_cp + phe$PC2_cp + 
                phe$PC3_cp + phe$PC4_cp + phe$PC5_cp + phe$PC6_cp + phe$PC7_cp + 
                phe$PC8_cp + phe$PC9_cp + phe$PC10_cp + phe$PC11_cp + phe$PC12_cp + 
                phe$PC13_cp + phe$PC14_cp + phe$PC15_cp + phe$PC16_cp + phe$PC17_cp + 
                phe$PC18_cp + phe$PC19_cp + phe$PC20_cp + phe$PC21_cp + phe$PC22_cp + 
                phe$PC23_cp + phe$PC24_cp + phe$PC25_cp + phe$PC26_cp + phe$PC27_cp + 
                phe$PC28_cp + phe$PC29_cp + phe$PC30_cp')
## PC1_cp are principal component of the control probes only necessary if you have used CPACOR to normalize your data.
# regression
phe=read.table("data/phe.txt",sep="\t",header=T) 
samples=as.character(phe$meth_ID)
samples=paste("X",samples,sep="")  # might not be necessary
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










