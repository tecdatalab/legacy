#!/usr/bin/perl

# Given a list of files provided thru the command line arguments
# It invokes the physics scoring program and outputs the score for each of those files

# Additionally, it is possible to provide a native structure to compare to and a flag that
# states that there is $filename.h that holds a version of the pdb with hydrogens added

if($#ARGV + 1 == 0)
{
	print "Usage: allscores.pl [--hydrogens] [--native native.pdb] pdb1 pdb2 ...\n";
}
else
{
	$paramindex = 0;
	$usehydrogens = 0;
	$native = "";
	if($ARGV[$paramindex] eq "--hydrogens")
	{
		$usehydrogens = 1;
		$paramindex++;
	}
	if($ARGV[$paramindex] eq "--native")
	{
		$native = $ARGV[$paramindex + 1];
		$paramindex += 2;
	}

	$first = 1;
	for(; $paramindex <= $#ARGV; $paramindex++)
	{
		$filename = $ARGV[$paramindex];
		$basename = `basename $filename`;
		chomp($basename);
		($pdbname) = split(/\./, $basename);
		$score_cmd = "score --decoy $filename";

		# add the hydrogens param or native if required
		if($usehydrogens)
		{
			$score_cmd .= " -h $filename.h";
		}
		if($native ne "")
		{
			# if no hydrogens file was indicated then use the same one
			if(!$use_hydrogens)
			{
				$score_cmd .= " -h $filename";
			}
			$score_cmd .= " --native $native";
		}

		@scorelines = `$score_cmd`;
		if($first)
		{
			$first = 0;
			print "PDB,";
			print $scorelines[0];
		}
		print "$pdbname,";
		print $scorelines[1];
	}
}
