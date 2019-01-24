#!/usr/bin/perl

# This script is intended to be used prior to the execution of Multi-LZerD, after the pairwise predictions have been generated
# It will keep only the top N predictions (with the highest scores). Note that the raw output from LZerD does not sort
# the output file, so this script will do that and keep the top predictions

if($#ARGV + 1 != 3)
{
	print "Assuming pairwise predictions have already been generated using LZerD, and the files are named following the <chain>-<chain>.out ";
	print "format, this script will sort the results and keep the top N predictions, overwriting the previous .out files.\n";
	print "The script receives 3 parameters: the directory where the output files are stored, the chain ID's used to create the pairwise ";
	print "predictions and the number of predictions that will be kept\n";
	print "\nUsage: multilzerd_sort_decoys <directory> <comma separated chain names> <top N>\n";
	print "Example multilzerd_sort_decoys ./predictions/1Z5S A,B,C,D 54000\n";
}
else
{
	# we need to create an intermediate sorted file; this will be its suffix
	$tmpfilesuffix = ".tmp";

	# parameters passed from the cmd line
	$directory = $ARGV[0];
	@chains = split(/,/, $ARGV[1]);
	$n = $ARGV[2];

	# there should be one output file for each combination of chain ID's
	foreach $rec_chain_index (0 .. $#chains)
	{
		foreach $lig_chain_index ($rec_chain_index + 1 .. $#chains)
		{
		
			$rec_chain = $chains[$rec_chain_index];
			$lig_chain = $chains[$lig_chain_index];

			$predictionfile = "$directory/$rec_chain-$lig_chain.out";
			$sortedfile = "$predictionfile$tmpfilesuffix";

			print "Sorting $rec_chain-$lig_chain.out...\n";

			# sort the file and keep the top N
			`sort -k 13 -n -r $predictionfile > $sortedfile`;
			`head -$n $sortedfile > $predictionfile`;

			# remove the tmp file
			`rm -f $sortedfile`;
		}
	}
	print "Finished sorting\n";
}
