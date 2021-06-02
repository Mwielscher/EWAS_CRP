#!/usr/bin/bash
#PBS -l select=1:ncpus=1:mem=5GB
#PBS -l walltime=03:00:00
#PBS -N run_chrSTATE

tid=$PBS_ARRAY_INDEX

###  15 jobs 

mkdir data
cp /rds/general/user/mwielsch/home/WORK/CRP/enrichment/ALL_loci/enrich_lists_ALL_MARKER_SD_distr.txt data/all_loci_SD_distribution.txt
cp /rds/general/user/mwielsch/home/WORK/CRP/enrichment/ALL_loci/enrich_lists_ALL_MARKER_loci.txt data/all_loci_significant.txt

cp /rds/general/user/mwielsch/home/WORK/CRP/enrichment/test_run/450K_SD_bins_for_permutation.RData data/
cp /rds/general/user/mwielsch/home/WORK/CRP/enrichment/data/anno.RData data/ 

ls data

module load R
R --vanilla --args ${tid} <<"EOF"

args=commandArgs(trailingOnly = TRUE)
tid=as.numeric(as.character(args[1]))
load("data/anno.RData")
hits=read.table("data/all_loci_significant.txt",sep=" ",header=T) 
#####   background CpGs
load("data/450K_SD_bins_for_permutation.RData")
input_dist=read.table("data/all_loci_SD_distribution.txt",header=T)
distr=input_dist$Freq

res=as.data.frame(matrix(NA,ncol=7,nrow=1))
colnames(res)=c("expected","expected_PERC","observed","observed_PERC","ePval_incr","ePval_decr","fPval")

all_trait=c("10_cell_E001","11_cell_E001","12_cell_E001","13_cell_E001","14_cell_E001","1_cell_E001","2_cell_E001","3_cell_E001","4_cell_E001","5_cell_E001","6_cell_E001","7_cell_E001","8_cell_E001","9_cell_E001")

trait=all_trait[tid]
rownames(res)=trait
ref.dat=read.table(paste("/rds/general/user/mwielsch/home/WORK/CRP/enrichment/data/FINAL/CHROMstateE001/","chromstate_",trait,".txt",sep=""),sep="\t")
ref=as.character(ref.dat$V1)

  obs=as.numeric(table(as.character(hits$MarkerName) %in% ref)[2])
	obs[is.na(obs)]=0 
 for (i in 1:10000){
	r1=as.character(bin1$V1[sample(1:nrow(bin1),distr[1])])
	r2=as.character(bin2$V1[sample(1:nrow(bin2),distr[2])])
	r3=as.character(bin3$V1[sample(1:nrow(bin3),distr[3])])
	r4=as.character(bin4$V1[sample(1:nrow(bin4),distr[4])])
	r5=as.character(bin5$V1[sample(1:nrow(bin5),distr[5])])
	r6=as.character(bin6$V1[sample(1:nrow(bin6),distr[6])])
	r7=as.character(bin7$V1[sample(1:nrow(bin7),distr[7])])
	r8=as.character(bin8$V1[sample(1:nrow(bin8),distr[8])])
	r9=as.character(bin9$V1[sample(1:nrow(bin9),distr[9])])
	r10=as.character(bin10$V1[sample(1:nrow(bin10),distr[10])])
	r11=as.character(bin11$V1[sample(1:nrow(bin11),distr[11])])
	r12=as.character(bin12$V1[sample(1:nrow(bin12),distr[12])])
	r13=as.character(bin13$V1[sample(1:nrow(bin13),distr[13])])
	rand=c(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13)

	perm=table(rand %in% ref)[2]
	perm[is.na(perm)]=0 
	if(exists("perm.dist")){
       perm.dist=c(perm.dist,perm)
     } else {
        perm.dist=perm
       }
    }
## calculate empirical p-value  
  p.vals_incr = sum(perm.dist >= obs) /10000
  p.val_decr = sum(perm.dist <= obs) /10000
#make Fisher test:
   all.cpg=c(416403)
  in.df=matrix(c(all.cpg-floor(mean(perm.dist)),all.cpg-obs,floor(mean(perm.dist)),obs),nrow=2)
  ## cont.table for Fisher test:
  ## all CpG - expected     expected
  ## all CpG - observed     observed
   fisher_test=tryCatch({fisher.test(in.df)$p.value},error = function(error) {return(rep(NA,1))})
  res[trait,]=c(mean(perm.dist),mean(perm.dist)/nrow(hits)*100,obs,obs/nrow(hits)*100,p.vals_incr,p.val_decr,fisher_test)
#  png(paste("/project/eph/matthias/CRP/analysis/annot_hits_2/result/plots/CHROM_state_",trait,".png",sep=""))
#  plot(density(perm.dist),xlim=c(min(perm.dist)-50,max(perm.dist)+50),main=paste("CRP-EWAS tophits enriched\n in ",trait,sep=""))
#  rect(obs,0,obs+5,0.01,col="red",border=NA)
#  dev.off()
  rm(perm.dist,fisher_test,p.vals_incr,p.val_decr,obs)


write.table(res,file=paste("/rds/general/user/mwielsch/home/WORK/CRP/enrichment/ALL_loci/result/res_CHROM_state_E001_",trait,".txt",sep=""),sep="\t",col.names=T,row.names=T,quote=F)


EOF

