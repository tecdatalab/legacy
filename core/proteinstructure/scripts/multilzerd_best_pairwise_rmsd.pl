#!/usr/bin/perl

# Taking info from the pairwise .stats files, this script will output the best LRMSD conformation
# for each of the files

if($#ARGV + 1 != 2)
{
	print "Taking info from the pairwise .stats files, this script will output the best LRMSD conformation";
	print "for each of the files\n";
	print "\nUsage: multilzerd_best_pairwise_rmsd.pl <directory> <comma separated chain names>\n\n";
	print "directory: Where the LZerD output files are located\n";
	print "chain names: Used to infer the names of LZerD output files in the directory\n";
	print "\nExample multilzerd_best_pairwise_rmsd.pl ./predictions/ A,B,C,D\n";
}
else
{
	$directory = $ARGV[0];
	@chains = split(/,/, $ARGV[1]);

	# generate a stats file for each pairwise prediction
	foreach $rec_chain_index (0 .. $#chains)
	{
		foreach $lig_chain_index ($rec_chain_index + 1 .. $#chains)
		{
		
			$rec_chain = $chains[$rec_chain_index];
			$lig_chain = $chains[$lig_chain_index];

			$stats_file = "$directory/$rec_chain-$lig_chain.stats";
			print "$rec_chain-$lig_chain\t";
			$best_rmsd =  `sort -k 3,3 -n $stats_file | head -1 | sed s/nan/0/g | awk '{print \$1 "\t" \$2 "\t" \$3 "\t" \$4}'`;
			chomp($best_rmsd);
			print "$best_rmsd\n";
		}
	}
}
