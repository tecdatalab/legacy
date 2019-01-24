#!/usr/bin/perl

# After calling "splitter", it creates a separate file for each chain/model
# Thus a model with 3 chains will create 3 separate files. If we want to
# integrate files to have 1 pdb per model we can call this script

if($#ARGV + 1 != 3)
{
	# start and end just determine the first index and the last index of the files to be processed
	print "Usage: joinchains.pl <split PDB's directory> <start> <end>\n";
}
else
{
	$directory = $ARGV[0];
	$start = $ARGV[1];
	$end = $ARGV[2];

	while($start <= $end)
	{
		$filenumber = sprintf("%04d", $start);
		`cat $directory/?-$filenumber*.pdb > $filenumber.pdb`;
		$start++;
	}
}
