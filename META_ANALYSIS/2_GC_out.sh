#!/bin/bash
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=01:20:00
#PBS -N GC_out



module load anaconda3/personal

R --vanilla << "EOF"

tt1=read.table("/project/eph/matthias/CRP/analysis/final_meta/COMPARE_IVW_JAN_allANCEST_het_check_11.tbl",sep="\t",header=T)
tt1=tt1[,!colnames(tt1) %in% c("Allele1","Allele2")]
head(tt1)
names(tt1)=c("MarkerName","Effect","StdErr","P_value", "Direction","HetISq","HetChiSq","HetDf","HetPVal")
tt2=read.table("/project/eph/matthias/CRP/analysis/final_meta/COMPARE_sampMETA_JAN_allANCEST_11.tbl",sep="\t",header=T)
tt3=tt2[,colnames(tt2) %in% c("MarkerName","Weight")]
head(tt3)
tt=merge(tt1,tt3,all.x=T,by.x="MarkerName",by.y="MarkerName")
head(tt)
 lambda = qchisq(median(tt$P_value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)	
 print(lambda)  
 tt$adj.P=pchisq((qchisq(tt$P_value, df=1, lower=F)/lambda), df=1, lower.tail=F) 

 
write.table(tt,file="/project/eph/matthias/CRP/analysis/final_meta/result/IVW_GC_OUT_HET_CHECK_2.txt",sep="\t",col.names=T,row.names=F,quote=F)

print(head(tt))



EOF
