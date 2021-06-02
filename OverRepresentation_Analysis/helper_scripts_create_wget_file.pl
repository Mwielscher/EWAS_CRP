#!/usr/bin/perl


use strict;
use warnings;
use WWW::Mechanize;
use Data::Dumper;

my $url = 'http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/gappedPeak/';
my $trait='H3K27me3.gappedPeak';


sub catch_them_all {
	my $mech  = WWW::Mechanize->new();
	$mech->get( $url );
	my @links = $mech->links();

	my @allFiles;
	foreach my $link (@links) {
		#print "LINK: ". $link->url() . "\n";
		push @allFiles, $link->url(); 
      	}
	print Dumper @allFiles;
	open (OUTPUT, ">","wget_roadmap.sh");
	print OUTPUT "#!/usr/bin/bash \n";
	foreach my $file(@allFiles) {
		print "loop = $file \n";
		if($file =~ /$trait/) 
		{
			my $print = join"","wget ",$url,$file,"\n";
			print OUTPUT $print;
		}
	}	
	close(OUTPUT);
	


}
catch_them_all();











