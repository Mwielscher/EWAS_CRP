### input data for this are produced by the script 1_assoc_CpG_cause_CRP.sh
rm(list=ls())
cpg=read.table("CpG_cause_CRP_instuments_INCL_DIRECT.txt",header=T,fill=T)
# we tested 8930738 sites i.e.5.598641e-09
# makes 8930738 genome wide test
# threshold is 5.6e-09 + 5e-08 makes 1x10-18  == 1e-18

cpg=cpg[cpg$P.value<1e-18 ,]
# direct effect cut in half
cpg$directHALFED=(-log10(cpg$P.value)-(-log10(cpg$DIRECT_P_val))) > ((-log10(cpg$P.value))/2)


## ------------------------------------------------------------------------------------------
#### TRIANGUATION
###--observed
## where beta3 == observed --> SNP ~ CRP
## beta 1 == SNP ~ CpG
## beta 2 == CRP ~ CpG   == EWAS result

cpg_tri=cpg[,colnames(cpg) %in% c("CpG_ID","Effect","StdErr","Effect2","StdErr.1")]
colnames(cpg_tri)=c("ID","beta1","beta1_SE" ,"beta3","beta3_SE")  #pred 1 = DNAmeth ~ SNP
cpg_tri$beta3=cpg_tri$beta3*-1
##--pred
meth1=meth[,colnames(meth) %in% c("MarkerName", "TE_Effect","TE_StdErr")]

cpg_tri=merge(cpg_tri,meth1,by.x="ID",by.y="MarkerName",all.x=T)
colnames(cpg_tri)[6:7]=c("beta2","beta2_SE")

cpg_tri$pred=cpg_tri$beta1*cpg_tri$beta2
cpg_tri$pred_SE = sqrt(cpg_tri$beta1_SE^2*cpg_tri$beta2_SE^2 + cpg_tri$beta1_SE^2 * cpg_tri$beta2^2 + cpg_tri$beta2_SE^2*cpg_tri$beta1^2)*10
cpg_tri$Z=cpg_tri$pred/cpg_tri$pred_SE
cpg_tri$cpg_Pval=2*pnorm(-abs(cpg_tri$Z))
plot(cpg_tri$pred, cpg_tri$beta3,xlim=c(-0.2,0.2),ylim=c(-0.2,0.2))

#---------------------------------------------------------------------------------------------------------
###########  ------------     make plot
cpg_tri=dat
#####  
plot(cpg_tri$pred[cpg_tri$ID %in% cpg_mr_hits],cpg_tri$observed[cpg_tri$ID %in% cpg_mr_hits])
summary(lm(cpg_tri$pred[cpg_tri$ID %in% cpg_mr_hits]~cpg_tri$observed[cpg_tri$ID %in% cpg_mr_hits]))
##############
xl=max(sqrt((cpg_tri$pred)^2))
yl=max(sqrt(cpg_tri$beta3^2))
xl=max(c(xl,yl))
yl=xl

plot(cpg_tri$pred, cpg_tri$beta3, xlim =c(-xl,xl), ylim=c(-yl,yl), xlab ='Predicted Effect', ylab='Observed Effect', col='black',bg='grey', lwd=2,pch=21, cex=1, main= "MR analysis: \n DNA methylation changes CRP levels", fg='grey', yaxt='n', xaxt='n')
axis(1,sort(unique(c(seq(0,-xl,-xl/3),seq(0,xl,xl/3)))), format(c((seq(-xl,-xl/3, xl/3)), 0.00, (seq(xl/3,xl, xl/3))),digits=1, nsmall=1), fg='grey', lwd.ticks=4)
axis(2,sort(unique(c(seq(0,-xl,-xl/3),seq(0,xl,xl/3)))), format(c((seq(-xl,-xl/3, xl/3)), 0.00, (seq(xl/3,xl, xl/3))),digits=1, nsmall=1), fg='grey',lwd.ticks=4, padj=0.5)
#abline(h=0, col='grey', lty=1, lwd=4)
abline(v=0, col='grey', lty=1, lwd=2)
abline(h=0, col='grey', lty=1, lwd=2)
abline(0,1, col= 'grey', lty=2, lwd=2)
points(cpg_tri$pred, cpg_tri$beta3, col='black',bg='blue', lwd=2,pch=21, cex=1)
############## significant
points(cpg_tri$pred[cpg_tri$color %in% c(2)], cpg_tri$beta3[cpg_tri$color %in% c(2)], col='black',bg='red', lwd=2,pch=21, cex=1)
box(fg='grey', lwd=3)



#---------------------------------------------------------------------------------------------------------
###  Mendelian Randomization: DNA meth   causes CRP changes---- ratio method --------
library(MendelianRandomization)
res2=as.data.frame(matrix(NA,nrow=(length(cpg$MarkerName)),ncol=4))
colnames(res2)=c("ID","IVW_Pvalue","IVW_Estimate","IVW_StdError")
i=0
for (k in as.character(cpg$CpG_ID)) { 
  print(c(k))
  i=i+1
  MR_cpg_cause = mr_input(bx = cpg$Effect[cpg$CpG_ID %in% c(k)],
                          by = cpg$Effect2[ cpg$CpG_ID %in% c(k)],
                          byse = cpg$StdErr.1[cpg$CpG_ID %in% c(k)])
  IVW= mr_ivw(MR_cpg_cause)
  res2[i,1:4]=c(as.character(k),IVW$Pvalue,IVW$Estimate,IVW$StdError)
  
}


