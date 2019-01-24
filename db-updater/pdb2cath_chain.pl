#!/usr/bin/perl

# Script: pdb2cath_chain.pl
# Description:
# Author: YI XIONG
# email: xiongyi128@gmail.com
# Date: Nov.11.18

my %cathhash;

$infile = "./pdbidlistupdate/cath-domain-list.txt";
open(CATH,$infile);
while (<CATH>) {
	next if /^#/;
	chomp($_);
	@data = split (/\s+/,$_);
	$id = substr($data[0],0,4)."-".substr($data[0],4,1);
	$cathhash{$id} .= "$data[1].$data[2].$data[3].$data[4],";
}
close(CATH);

$infile = "./pdbidlistupdate/chainidlist.txt";
open(IN,$infile);
@idlistarray=<IN>;
close(IN);
chomp(@idlistarray);

$outfile = "./pdbidlistupdate/chainid.cath.mapping.txt";
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
			print OUT $key.",";	
		}
		print OUT "\n";
	}
}
close(OUT);

print "Success\n";

