#!/usr/bin/perl

# Assuming we have a multiple docking prediction file and the original PDB, this script performs 2 processes
# 1. Re-create the PDB file for each of the predictions and calculates the Zernike Descriptors that define its shape
# 2. Using those PDB files, it creates a simulated EM representation using the EMAN2 package and also generates
# 	the Zernike Descriptors that define their shape.
# (ZD's are generated for the native PDB and the EM simulation of the native PDB too)
# The purpose of this script is to generate data that allows us to determine if 3DZD's can be used to compare the shape
# of EM images.

# this procedure converts directly from a PDB to inv representation (without generatinf the EM image)
sub pdb2inv
{
	my ($prefix, $directory) = @_;
	`zernike.pl $directory$prefix.pdb`;
}

# provide the prefix of a pdb file in order to generate a simualted EM image first, and then
# the inv file with the descriptors corresponding to that EM image
# the last parameter determines the name given to the output file
sub pdb2em2inv
{
	my ($prefix, $low, $high, $dir) = @_;

	my $pdb_file = "$dir$prefix.pdb";
	my $mrc_file = "$prefix.mrc";
	my $situs_file = "$prefix.situs";
	my $grid_file = "$prefix.grid";
	my $inv_file = "$prefix.inv";
	`e2pdb2mrc.py --apix=0.5 --res=10 $pdb_file $mrc_file`;
	#`~/bin/EMAN2/bin/e2pdb2mrc.py $pdb_file $mrc_file`;
	# map2map is interactive through the command line. Option 2 specifies the program to use MRC files
	# as input. Since it receives it as standard input then a "2" is piped as std input to the command
	`echo 2 | map2map $mrc_file $situs_file`;
	`situs2grid.pl $situs_file $grid_file $low $high`;
	`zernike $grid_file 20`;
	
	#delete all files except the inv
	`rm $mrc_file $situs_file $grid_file`;
	`mv $inv_file $prefix.em.inv`;
}

if($#ARGV + 1 != 5)
{
	print "Usage: mdga2inv.pl <PDBID> <number of predictions> <em lower bound> <em upper bound> <mdga input dir>\n\n";
	print "The script assumes that there is a <PDBID>.ga.out file and a <PDBID>.pdb in the current directory\n";
}
else
{
	$pdbid = $ARGV[0];
	$number_of_predictions = $ARGV[1];
	$lower_bound = $ARGV[2];
	$upper_bound = $ARGV[3];
	$input_dir = $ARGV[4];

	$euclidean_distance_file = "$pdbid.euc";

	# generate the zernike descriptors file (INV) for both EM and PDB representations
	pdb2em2inv($pdbid, $lower_bound, $upper_bound, "./");
	pdb2inv($pdbid, "./");

	# the first line of the euclidean distances file is the distance between these 2 inv files
	$nativepdb2nativeem = `invdistance.pl $pdbid.em.inv $pdbid.pdb.inv`;
	chomp($nativepdb2nativeem); 
	open(EUC_FILE, ">$euclidean_distance_file");
	print(EUC_FILE "Native PDB to Native EM distance: $nativepdb2nativeem\n");
	print(EUC_FILE "\nDecoy: PDB-PDB distance\tEM-EM distance\n");
	close(EUC_FILE);

	# and now proceed to the decoys...
	
	# create the pdb files from the .ga.out file
	$mdga_file = "$pdbid.ga.out";
	`mdga_create_pdb $mdga_file $input_dir 1 $number_of_predictions $pdbid`;


	# 2. for each file generate the descriptors file
	foreach $decoy_file(`ls $pdbid-*.pdb`)
	{
		chomp($decoy_file);
		@split_filename = split(/\.pdb/, $decoy_file);
		$prefix = $split_filename[0];
		# generate zernike descriptors for both PDB and EM representations

		# "./" because, for now, they are generated in the current dir
		pdb2em2inv($prefix, $lower_bound, $upper_bound, "./");
		pdb2inv($prefix, "./");
		# Calculate euclidean distance from
		# 1. Native pdb.inv to Decoy pdb.inv
		$pdbdistance = `invdistance.pl $pdbid.pdb.inv $prefix.pdb.inv`;
		chomp($pdbdistance);
		# 2. Native em.inv to Decoy em.inv
		$emdistance = `invdistance.pl $pdbid.em.inv $prefix.em.inv`;
		chomp($emdistance);
		# Write the distances to the .euc file
		# The file is opened and closed each time because each iteration of the loop can take a considerable time
		open(EUC_FILE, ">>$euclidean_distance_file");
		print(EUC_FILE "$prefix\t$pdbdistance\t$emdistance\n");
		close(EUC_FILE);
	}
}
