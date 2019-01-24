#!/usr/bin/perl

use strict;

my $pdb_file = $ARGV[0];
my $type = $ARGV[1];
my ($pdb_id) = $pdb_file =~ /(\w+)\.pdb/;

my $bin = './zd';
my $msroll = "$bin/msroll";
my $poly2obj = "$bin/poly2obj.pl";
my $msms = "$bin/msms";
my $pdb2xyzr = "$bin/pdb_to_xyzr";
my $msms2any = "$bin/msms2any";
my $vox = "$bin/obj2grid";
my $zernike = "$bin/map2zernike";

my $time_start = time();
my $run_time;


if ($type eq 'msms') {

	# Covert PDB to XYZR format
	print "$pdb2xyzr $pdb_file > $pdb_file.xyzr\n";
	`$pdb2xyzr $pdb_file > $pdb_file.xyzr`;
	
	# Triangulate Molecular Surface by MSMS
	print "---- Surface Triangulation by MSMS ----\n";
	print "$msms -if $pdb_file.xyzr -of $pdb_file\n";
	`$msms -if $pdb_file.xyzr -of $pdb_file`;
	$run_time = time() - $time_start;
	print "Run Time: $run_time sec\n";

	# Convert MSMS Format to OBJ Format
	`$msms2any $pdb_file.face $pdb_file.vert 0`;

} 

unless (-s "$pdb_file.vert") {

	# Triangulate Molecular Surface by MSP
	print "---- Surface Triangulation by MSP ----\n";
	`$msroll -m $pdb_file -t $pdb_file.poly`;
	$run_time = time() - $time_start;
	print "Run Time: $run_time sec\n";

	# Convert MSP Poly Format to OBJ Format
	`$poly2obj $pdb_file.poly > $pdb_file.obj`;
}


# Voxelize
my $vox_time_start = time();
print "---- Surface Voxelization ----\n";
`$vox -g 64 $pdb_file.obj`;
my $vox_run_time = time() - $vox_time_start;
print "Run Time: $vox_run_time sec\n";

# 3D Zernike Transform
my $zernike_time_start = time();
print "---- 3D Zernike Transform ----\n";
`$zernike $pdb_file.obj.grid -c 0.5`;
my $zernike_run_time = time() - $zernike_time_start;
`mv $pdb_file.obj.grid.inv $pdb_file.inv`;
print "Output INV: $pdb_id\.inv\n";
print "Run Time: $zernike_run_time sec\n";

# Clean up
if ($type eq 'msms') {
	`rm $pdb_file.xyzr`;
	`rm $pdb_file.face`;
	`rm $pdb_file.vert`;
} 

if (-s "$pdb_file.poly") {
	`rm $pdb_file.poly`;
}

`rm $pdb_file.obj`;
`rm $pdb_file.obj.grid`;

print "\n---- Completed ----\n";
my $total_run_time = time() - $time_start;
print "\nTotal Run Time: $total_run_time sec\n\n";

