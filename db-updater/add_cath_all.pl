#!/usr/bin/perl

# Script: add_cath.pl
# Description:
# Author: YI XIONG
# email: xiongyi128@gmail.com
# Date: Nov.11.19

$pdb2cath        = $ARGV[0];
$input_filename  = $ARGV[1];
$output_filename = $ARGV[2];

%id_cath_hash;
open( IN, $pdb2cath );
while (<IN>) {

	chomp($_);
	@data = split( /\t/, $_ );

	$pdb_id                = $data[0];
	@catharray             = split( /,/, $data[1] );
	$id_cath_hash{$pdb_id} = $catharray[0];
}
close(IN);

# read from the basic pdb_latest file in order to add to it
open( IN_BASIC, $input_filename );

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
	if ( exists $id_cath_hash{$pdb_id} ) {

		$cath = $id_cath_hash{$pdb_id};
		( $c, $a, $t, $h ) = $cath =~ /([0-9]*)\.([0-9]*)\.([0-9]*)\.([0-9]*)/;
		print NEW_PDB_LATEST "$line $c $a $t $h\n";
	}
	else {
		print NEW_PDB_LATEST "$line -1 -1 -1 -1\n";
	}
}

#close both files
close(NEW_PDB_LATEST);
close(IN_BASIC);
