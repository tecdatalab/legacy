#!/usr/bin/perl

# Script: pdb2cath_complex.pl
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
	$id = substr($data[0],0,4);
	$cathhash{$id} .= "$data[1].$data[2].$data[3].$data[4],";
}
close(CATH);

$infile = "./pdbidlistupdate/complexidlist.txt";
open(IN,$infile);
@idlistarray=<IN>;
close(IN);
chomp(@idlistarray);

$outfile = "./pdbidlistupdate/complexid.cath.mapping.txt";
open(OUT,">$outfile"); 
for my $i (0..$#idlistarray) {
	$pdb_id = substr($idlistarray[$i], 0, 4);
	if(exists $cathhash{$pdb_id}){
		my %tmplisthash;
		my @tmparray=split (/,/, $cathhash{$pdb_id});
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

