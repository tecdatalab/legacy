#!/usr/bin/perl

# Usage: poly2obj.pl <msp_output>
# This program converts MSP POLY output to OBJ format
# Written by David La (2006)

use strict;

my @line;
my $vert_count;
my $face_count;
my $col_count;

print "#\r\n#\r\n";
print "mtllib $ARGV[0]\r\n";
print "# object BLACK TRIM\r\n";
print "g BLACK_TRIM\r\n";

while (<>) {
	
	chomp;
	s/^\s+//;
	
	@line = split /\s+/;

	$col_count = scalar @line;

	if ( $col_count == 12) {
		$vert_count++;
		print "v $line[0] $line[1] $line[2]\r\n";
	}
	elsif ($col_count == 9) {
		$face_count++;
		if ($face_count == 1) {
			print "# $vert_count vertices\r\n";
			print "usemtl BLACK_TRIM\r\n";
		}
		print "f  $line[3] $line[4] $line[5]\r\n";
	}
}


