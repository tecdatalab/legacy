#!/usr/bin/perl

# Script: pdb2cath_domain.pl
# Description:
# Author: YI XIONG
# email: xiongyi128@gmail.com
# Date: Nov.11.13
# edit 6.29.2015: not include domains that have been removed

my %cathhash;

$infile = "./pdbidlistupdate/cath-domain-list.txt";
open(CATH,$infile);
while (<CATH>) {
	next if /^#/;
	chomp($_);
	@data = split (/\s+/,$_);
	$id = substr($data[0],0,4)."-".substr($data[0],4,1)."-".substr($data[0],5,2);

	my $two = substr($data[0],5,2);
	if ($two ne "00") {
		$cathhash{$id} .= "$data[1].$data[2].$data[3].$data[4],";
	}
}
close(CATH);

$infile = "./pdbidlistupdate/domainidlist.txt";
open(IN,$infile);
@idlistarray=<IN>;
close(IN);
chomp(@idlistarray);

$outfile = "./pdbidlistupdate/domainid.cath.mapping.txt";
open(OUT,">$outfile"); 

for my $i (0..$#idlistarray) {
	if(exists $cathhash{$idlistarray[$i]}){
		my %tmplisthash;
		my @tmparray=split (/,/, $cathhash{$idlistarray[$i]});
		for my $j (0..$#tmparray) {
			$tmplisthash{$tmparray[$j]} = 1;	
		}
		
		print OUT $idlistarray[$i]."\t";
		foreach $key ( sort keys %tmplisthash ) {
			print OUT $key;	
		}
		print OUT "\n";
	}
}
close(OUT);

print "Success\n";

