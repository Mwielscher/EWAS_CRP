#!/usr/bin/perl

use strict;
use warnings;

my $dir="/project/eph/matthias/CRP/annot_hits_2/data/HiC/input_data/";
my @traits=("A1","A2","B1","B2","B3","B4");



sub make_list {
	#my$trait='GSE63525_HeLa_HiCCUPS_looplist.txt';
	foreach my$trait(@traits) {
		open(IN,$dir."GSE63525_GM12878_subcompartments.bed") or die "cannot open $dir $trait IN";
		open(OUT,">",$dir."final_list_compartment_".$trait) or die "cannot open out";
		<IN>;
		while(<IN>){
			chomp $_;
			my@line=split "\t",$_;
				if($line[3] eq $trait) {
			print OUT "$line[0] $line[1] $line[2] \n";
			} else { next; }
		}
		close(IN);
		close(OUT);
	}
}
make_list();
