#!/usr/bin/perl


use strict;
use warnings;
use Data::Dumper;
my $dir='/project/eph/matthias/CRP/analysis/annot_hits_2/data/STATEBYLINE/';
my $sel='E128';
my @chrom= (1..22);
sub make_bins {
       	foreach my$chr(@chrom) {
		 open(IN,$dir."input_data/".$sel."_15_coreMarks_chr".$chr."_statebyline.txt") or die "cannot open IN";
		open(OUT,">",$dir.$sel."/15coreMarks_chr".$chr."_".$sel.".txt") or die "cannot open OUT";
		<IN>;
		my $i = 0;
		while(<IN>) {
		chomp;
		if ($_=~/^Max/) { next; }
			else {
			my $k= $i +200;
			my $s= $i +1;
			print OUT "chr$chr $s $k $_ \n";
			$i = $k;
			}
		}
		close(IN);
		close(OUT);
	}
}

make_bins;



