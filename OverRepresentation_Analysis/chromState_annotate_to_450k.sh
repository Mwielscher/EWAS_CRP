#!/bin/bash

dir='/project/eph/matthias/CRP/analysis/annot_hits_2/data/'
final='FINAL/'
traits=$(echo{1..14})
name='CORE_'
subdir='STATEBYLINE/E128/'
newname='CHROMstate'
sel='E128'
#################################################   Genome Build hg19
mkdir $dir$final$newname$sel
for trait in $traits
do
cat <<EOF >"runCORE_$trait.pl"
#!/usr/bin/perl
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -N CORE_$trait
use strict;
use warnings;
use Data::Dumper;

my \$dir = "$dir";
sub match_it {
	my @ref;
	open(REF,\$dir."450k_chromPOS_both.txt") or die "cannot open REF";
	<REF>;
	while (<REF>) {
		chomp;
		my@line = split "\t",\$_;	
		my\$vec= join " ",\$line[0],\$line[1],\$line[2];                 #### HG19 !!!!!!!!!
		push @ref,\$vec;
	}
	open(IN,"${dir}${subdir}${name}${trait}"."_${sel}.txt") or die "cannot open IN";
	open(OUT,">","${dir}${final}${newname}${sel}"."/"."chromstate_${trait}"."_cell_${sel}.txt") or die "cannot open OUT";

	<IN>;
	while(<IN>) {
		chomp;
		my @line=split " ",\$_;	
		my\$chr =\$line[0];
		my\$start=\$line[1];
		my\$end=\$line[2];
		foreach my\$in (@ref) {
			my@illum = split " ",\$in;
			if (\$illum[1] eq \$chr && \$illum[2] > \$start && \$illum[2] < \$end ) {
			print OUT "\$illum[0]\n";		
			}else {
			next; } 	
		}

	}
	close(IN);
	close(OUT);
}

match_it;
EOF

chmod +x runCORE_$trait.pl
sleep 0.01
qsub runCORE_$trait.pl
done



