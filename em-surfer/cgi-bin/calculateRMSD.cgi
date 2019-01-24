#!/usr/bin/perl

use CGI;

$q = new CGI;

my @params = $q->param();
# the 0th parameter is the pdb id of the given query
# 1-25th parameters are ids of comparison pdb
my $a1 = $q->param(@params[0]);
my $a2;
my $i;
my $ret = "rmsd";
my $ce = "";
my $pdb_chains_dir = "../../db/chain/";
my ($size, $alignlen, $coverage);

print $q->header;

for ($i = 1; $i <= 25; $i++) {
	
	$a2 = $q->param(@params[$i]);
	if ($a2 != 0) {
	
		$line = `cd ../ce_distr; ./CE - $pdb_chains_dir$a1.pdb - $pdb_chains_dir$a2.pdb - ./scratch > ../cgi-bin/tmp/ce\_$$\_$i.txt; cat ../cgi-bin/tmp/ce\_$$\_$i.txt`;
		$ce = $ce . "`;" . $line;

		($rmsd) = $line =~ /Rmsd = (\d+\.\d+)A/;
		$rmsd = 'N/A' if !$rmsd;
		$ret = $ret . "`;" . $rmsd;
		
		($size) = $line =~ /Size=(\d+)/;
		($alignlen) = $line =~ /Alignment length = (\d+)/;
		$coverage=$alignlen/$size*100;
		$coverage = " (".sprintf("%.1f", $coverage) . "%".")";

	} else {
		
		$ce = $ce . "`;0";
		$ret = $ret . "`;0";
	}
}

print "$ret$ce`;$a1`;$$`;$coverage";
