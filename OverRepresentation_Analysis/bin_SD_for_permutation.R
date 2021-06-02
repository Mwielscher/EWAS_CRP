##############   bin background list
back=read.table("prelim_binning_input.txt",header=T,fill=T)
colnames(back)[2]=c("SD")
back$SD=as.numeric(as.character(back$SD))
back$combMEAN=as.numeric(as.character(back$combMEAN))
back=back[!is.na(back$SD) ,]
back=back[!is.na(back$combMEAN) ,]
quantile(back$SD,probs=seq(0,1,0.05))
back$SD[back$SD<=0.0069]=0.007
back$SD[back$SD>=0.075]=0.075
back$catSD=cut(back$SD,seq(0.007,0.075,0.005),right=F,labels=c(1:13))

for (i in 1:13) {
write.table(as.data.frame(back$probeID[back$catSD %in% c(i)]),file=paste("binned_SD_all_",i,sep=""),col.names=F,row.names=F,quote=F)
}
## we also created an .RData containing bin1 to bin13 for the enrichment scripts

##### create SD bin distriubtions of signifcant loci
dat=read.table("Signif_loci_for_enrichment_gene_mean_SD.txt",header=T)  #
inter=back[,colnames(back) %in% c("probeID","catSD")]
dat1=merge(dat,inter,by.x="MarkerName",by.y="probeID",all.x=T)
dat1=dat1[!is.na(dat1$catSD),]
input_dist_all=as.data.frame(table(dat1$catSD))
write.table(input_dist_all,file="all_loci_SD_distribution.txt",col.names=T,row.names=F,quote=F)
write.table(dat1,file="all_loci_significant.txt",col.names=T,row.names=F,quote=F)





