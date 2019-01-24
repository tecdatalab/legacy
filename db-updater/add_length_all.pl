#!/usr/bin/perl

# Script: add_length.pl
# Description: Usage: add_length_chain.pl <chain directory> <basic pdb_latest file> <output filename>
# Author: YI XIONG
# email: xiongyi128@gmail.com
# Date: Nov.11.19

$pdb2length_file           = $ARGV[0];
$basic_pdb_latest_filename = $ARGV[1];
$output_filename           = $ARGV[2];

%id_length_hash;
open( IN, $pdb2length_file );
while (<IN>) {

	chomp($_);
	@data = split( /\t/, $_ );
	$id_length_hash{ $data[0] } = $data[1];
}
close(IN);

# read from the basic pdb_latest file in order to add to it
open( IN_BASIC, $basic_pdb_latest_filename );

open( NEW_PDB_LATEST, ">$output_filename" );

#read and re-write the header (2 lines)
$line = <IN_BASIC>;
print NEW_PDB_LATEST $line;
$line = <IN_BASIC>;
print NEW_PDB_LATEST $line;

while (<IN_BASIC>) {

	$line = $_;
	chomp($line);

	$pdb_id = substr( $line, 0, index( $line, ".pdb.inv" ) );
	if ( exists $id_length_hash{$pdb_id} ) {

		$length = $id_length_hash{$pdb_id};
		print NEW_PDB_LATEST "$line $length\n";
	}
	else {

		print NEW_PDB_LATEST "$line 0\n";
	}
}

#close both files
close(NEW_PDB_LATEST);
close(IN_BASIC);
