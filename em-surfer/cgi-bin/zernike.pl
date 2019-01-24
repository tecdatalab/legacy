#!/usr/bin/perl

use strict;

my $pdb_file = $ARGV[0];
# the directory where the executables and scripts are must be passed as parameter
my $bin = $ARGV[1];
my ($pdb_id) = $pdb_file =~ /(\w+)\.pdb/;

my $msroll = "$bin/msroll";
my $poly2obj = "$bin/poly2obj.pl";
my $msms = "$bin/msms";
my $pdb2xyzr = "$bin/pdb_to_xyzr";
my $msms2any = "$bin/msms2any";
my $vox = "$bin/obj2grid";
my $zernike = "$bin/zernike";

my $time_start = time();
my $run_time;

# Triangulate Molecular Surface by MSP
print "---- Surface Triangulation by MSP ----\n";
print `$msroll -m $pdb_file -t $pdb_file.poly`;

if (-s "$pdb_file.poly") {
	$run_time = time() - $time_start;
	print "Run Time: $run_time sec\n";
} else {
	print "Run Time: Did not complete ... running MSMS instead.\n";
}

# Convert MSP Poly Format to OBJ Format
print `$poly2obj $pdb_file.poly > $pdb_file.obj`;

# If MSP fails, then use MSMS
unless (-s "$pdb_file.poly") {
	# Covert PDB to XYZR format
	print "$pdb2xyzr $pdb_file > $pdb_file.xyzr\n";
	print `$pdb2xyzr $pdb_file > $pdb_file.xyzr`;
	
	# Triangulate Molecular Surface by MSMS
	print "---- Surface Triangulation by MSMS ----\n";
	print "$msms -if $pdb_file.xyzr -of $pdb_file\n";
	print `$msms -if $pdb_file.xyzr -of $pdb_file`;
	$run_time = time() - $time_start;
	print "Run Time: $run_time sec\n";

	# Convert MSMS Format to OBJ Format
	print `$msms2any $pdb_file.face $pdb_file.vert 0`;

}

# Voxelize
my $vox_time_start = time();
print "---- Surface Voxelization ----\n";
print `$vox -g 64 $pdb_file.obj`;
my $vox_run_time = time() - $vox_time_start;
print "Run Time: $vox_run_time sec\n";

# 3D Zernike Transform
my $zernike_time_start = time();
print "---- 3D Zernike Transform ----\n";
print `$zernike $pdb_file.obj.grid 20`;
my $zernike_run_time = time() - $zernike_time_start;
`mv $pdb_file.obj.inv $pdb_file.inv`;
print "Output INV: $pdb_id\.inv\n";
print "Run Time: $zernike_run_time sec\n";

# Clean up

`rm $pdb_file.poly`;
`rm $pdb_file.obj`;
`rm $pdb_file.obj.grid`;

if (-f "$pdb_file.xyzr") {
	`rm $pdb_file.xyzr`;
	`rm $pdb_file.face`;
	`rm $pdb_file.vert`;
}

print "\n---- Completed ----\n";
my $total_run_time = time() - $time_start;
print "\nTotal Run Time: $total_run_time sec\n\n";

