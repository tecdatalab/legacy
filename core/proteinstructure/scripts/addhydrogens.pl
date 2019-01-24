#!/usr/bin/perl
if($#ARGV + 1 < 1)
{
	print "Given a file.pdb generates a file.pdb.h with hydrogens added to it.\n";
	print "\nUsage: addhydrogens.pl <PDB file1> <PDB file2> ...\n";
}
else
{
	for $filename(@ARGV)
	{
		$basename = `basename $filename`;
		chomp($basename);
		@basename_parts = split(/\./, $basename);
		$prefix = $basename_parts[0];
	
		`hbplus -o $filename ; modify_hbplus $prefix.h $prefix.pdb.h ; rm $prefix.h $prefix.hb2`;
	}
}
