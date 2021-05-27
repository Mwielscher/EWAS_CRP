
## in script 1 we checked the effect sizes - then we performed a standard meta-analysis 
# with METAL software to create the input file for this analsis 
rm(list=ls())
library(qqman)
crp=read.table("CRP_causes_DNAmeth_META_analysed_FINAL_1.txt",header=T,fill=T)
#####  -------------- subset
cpg=read.table("CpG_cause_CRP_instuments_INCL_DIRECT.txt",header=T,fill=T)
# we tested 8930738 sites i.e.5.598641e-09
# makes 8930738 genome wide test
# threshold is 5.6e-09 + 5e-08 makes 1x10-18  == 1e-18
cpg=cpg[cpg$P.value<1e-18 ,]
###
crp=crp[crp$MarkerName %in% cpg$CpG_ID,]
##################  conventional MR

library(MendelianRandomization)
#install.packages("MendelianRandomization")
res2=as.data.frame(matrix(NA,nrow=(length(crp$MarkerName)),ncol=4))
colnames(res2)=c("ID","IVW_Pvalue","IVW_Estimate","IVW_StdError")
i=0
for (k in as.character(crp$MarkerName)) { 
  print(c(k))
  i=i+1
  MR_crp_cause = mr_input(bx = as.numeric(crp.assoc[1]),
                          by = crp$Effect[ crp$MarkerName %in% c(k)],
                          byse = crp$StdErr[crp$MarkerName %in% c(k)])
  IVW= mr_ivw(MR_crp_cause)
  res2[i,1:4]=c(as.character(k),IVW$Pvalue,IVW$Estimate,IVW$StdError)
  
}

write.csv(res2,file="CRP_causes_CpG_FIANL.csv")

###############   triangulation
#### TRIANGUATION
###--observed
## where beta3 == observed --> GRScore ~ CpG
## beta 1 == GRScore ~ CRP
## beta 2 == CRP ~ CpG   == EWAS result
crp_tri=crp[,colnames(crp) %in% c("MarkerName","Effect", "StdErr")]
colnames(crp_tri)=c("ID" ,"beta3","beta3_SE")  
crp_tri$beta1=rep(as.numeric(crp.assoc[1]),nrow(crp_tri))
crp_tri$beta1_SE=rep(as.numeric(crp.assoc[2]),nrow(crp_tri))

#cpg_tri$beta3=cpg_tri$beta3*-1
##--pred
meth1=meth[,colnames(meth) %in% c("ID", "TE_Effect","TE_StdErr")]

crp_tri=merge(crp_tri,meth1,by.x="ID",by.y="ID",all.x=T,all.y=T)
colnames(crp_tri)[6:7]=c("beta2","beta2_SE")
crp_tri=crp_tri[!is.na(crp_tri$beta3),]

crp_tri$pred=crp_tri$beta1*crp_tri$beta2
crp_tri$pred_SE = sqrt(crp_tri$beta1_SE^2*crp_tri$beta2_SE^2 + crp_tri$beta1_SE^2 * crp_tri$beta2^2 + crp_tri$beta2_SE^2*crp_tri$beta1^2)*10
crp_tri$Z=crp_tri$pred/crp_tri$pred_SE
crp_tri$crp_Pval=2*pnorm(-abs(crp_tri$Z))

cor.test(crp_tri$pred,crp_tri$beta3)   


###   ---------------------  create plot
xl=max(sqrt((crp_tri$pred/1000)^2))
yl=max(sqrt(crp_tri$beta3^2))
xl=max(c(xl,yl))
yl=xl

plot(crp_tri$pred/1000, crp_tri$beta3, xlim =c(-xl,xl), ylim=c(-yl,yl), xlab ='Predicted Effect', ylab='Observed Effect', col='black',bg='grey', lwd=2,pch=21, cex=1, main= "MR analysis: \n DNA methylation is consequential of CRP levels", fg='grey', yaxt='n', xaxt='n')
axis(1,sort(unique(c(seq(0,-xl,-xl/3),seq(0,xl,xl/3)))), format(c((seq(-xl,-xl/3, xl/3)), 0.00, (seq(xl/3,xl, xl/3))),digits=1, nsmall=1), fg='grey', lwd.ticks=4)
axis(2,sort(unique(c(seq(0,-xl,-xl/3),seq(0,xl,xl/3)))), format(c((seq(-xl,-xl/3, xl/3)), 0.00, (seq(xl/3,xl, xl/3))),digits=1, nsmall=1), fg='grey',lwd.ticks=4, padj=0.5)
#abline(h=0, col='grey', lty=1, lwd=4)
abline(v=0, col='grey', lty=1, lwd=2)
abline(h=0, col='grey', lty=1, lwd=2)
abline(0,1, col= 'grey', lty=2, lwd=2)
points(crp_tri$pred/1000, crp_tri$beta3, col='black',bg='blue', lwd=2,pch=21, cex=1)
box(fg='grey', lwd=3)
crp_tri$check=crp_tri$beta3*crp_tri$pred
table(crp_tri$check>0)
cor.test(crp_tri$pred/1000, crp_tri$beta3)
binom.test(table(crp_tri$check>0)[2], nrow(crp_tri)) 


write.csv(res2,file="MR_result_CRP_causes_DNA_meth.txt",quote=F)
write.csv(crp_tri,file="TRIANGULATION_CRP_causes_DNA_meth.txt",quote=F)
