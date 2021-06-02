#!/usr/bin/perl
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -N createTF

use strict;
use warnings;


my $dir='/project/eph/matthias/CRP/annot_hits_2/data/TF_encode/';
my @allTF=("ZBTB33","CEBPB","CTCF","TAF1","GABPA","USF1","SP1","EGR1","FOXA1","RUNX3","MAZ","RAD21","SMC3","MAFF","MAFK","BHLHE40","FOSL2","JUND","E2F6","MAX","POLR2A","PAX5","PHF8","PML","YY1","SIN3AK20","E2F1","GTF2F1","ATF2","MYC","KDM5A","MXI1","POU2F2","KDM5B","TBP","IRF1","EP300","TAF7","ELK1","RFX5","TCF7L2","CHD2","FOXP2","ATF3","BRCA1","NFYA","RELA","NFYB","GRp20","REST","JUN","E2F4","SRF","ELF1","CREB1","ATF1","SIX5","USF2","FOS","TBL1XR1","ZNF143","SP2","EBF1","CTCFL","TEAD4","THAP1","ZEB1","ZNF263","PBX3","UBTF","CBX3","BCLAF1","NR2C2","RBBP5","GATA1","RCOR1","FOSL1","GATA2","TAL1","GATA3","TCF12","BCL3","NFATC1","MEF2A","MEF2C","CCNT2","BACH1","HDAC2","TCF3","ZNF274","STAT1","BATF","SPI1","HMGN3","SETDB1","ETS1","ZBTB7A","EZH2","JUNB","SP4","TFAP2A","TFAP2C","NR2F2","ESR1","SIN3A","TRIM28","HNF4G","RXRA","GTF3C2","SUZ12","CTBP2","NR3C1","SAP30","CHD1","KAP1","NANOG","STAT5A","HDAC1","ELK4","NRF1","STAT3","HNF4A","FOXA2","SMARCC1","SMARCB1","ESRRA","STAT2","MYBL2","NFIC","SREBP1","ARID3A","CEBPD","IRF4","BCL11A","MTA3","FOXM1","ZNF217","HSF1","HDAC8","NFE2","IRF3","WRNIP1","GTF2B","HDAC6","SMARCA4","ZKSCAN1","BRF2","IKZF1","POU5F1","RPC155","PPARGC1A","BDP1","SIRT6","SMARCC2","MBD4","PRDM1","FAM48A","RDBP","ZZZ3","POLR3G","BRF1"
);

sub make_core {
		#open(IN,$dir."/input_data/wgEncodeRegTfbsClusteredV3.bed") or die "cannot open IN";
		foreach my$TF (@allTF) {
			print "doing $TF now !!\n";
			open(OUT,">",$dir."/input_data/".$TF."encode_hg19.txt") or die "cannot open OUT $TF";
			 open(IN,$dir."/input_data/wgEncodeRegTfbsClusteredV3.bed") or die "cannot open IN";
			while(<IN>) {
				chomp;
				my@line= split "\t",$_;
				if ($line[3] eq $TF ) {
				print OUT "$_\n";	
				}else { next; }
			}	
			close(OUT);
			close(IN);
		}
		
		
}
make_core();
