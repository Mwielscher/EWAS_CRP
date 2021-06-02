#!/usr/bin/perl
#PBS -l walltime=2:00:00
#PBS -lselect=1:ncpus=1:mem=10GB
#PBS -N CORE_

use strict;
use warnings;

my $sel='E128';
my $dir='/project/eph/matthias/CRP/analysis/annot_hits_2/data/STATEBYLINE/';
my @chrom= (1..22);
my @cores=(1..14);

sub make_core {
	foreach my$core(@cores) {
		open(OUT,">",$dir.$sel."/CORE_".$core."_".$sel.".txt") or die "cannot open OUT";
		foreach my$chr (@chrom) {
			open(IN,$dir.$sel."/15coreMarks_chr".$chr."_".$sel.".txt") or die "cannot open IN $chr";
			while(<IN>) {
				chomp;
				my@line= split " ",$_;
				if ($line[3] eq $core ) {
				print OUT "$_\n";
			#	print "$_ \n";
				}else { next; }
			}	
			close(IN);
		}
		close(OUT);
	}
}
make_core();
