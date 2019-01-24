#!/usr/bin/perl

# This script will trigger the clustering of pairwise predictions that have been generated already for
# a multi-chain complex

if($#ARGV + 1 != 8)
{
	print "Assuming pairwise LZerD predictions have been generated already, and that they follow the <chain>-<chain>.out naming format, ";
	print "this script will initiate the execution of the pairwise lzerd_clustering program for each pair of chains.\n";
	print "\nUsage: multilzerd_pairwise_cluster.pl <directory> <PDB ID> <comma separated chain names> <suffix> <RMSD cutoff> <allrmsd/lrmsd> ";
	print "<nosizes/sizes> <execute/print>\n\n";
	print "directory: Where the LZerD output files are located\n";
	print "PDB ID and chain names: Used to infer the names of PDB and LZerD output files in the directory\n";
	print "suffix: the output file will be of the form <chain>-<chain>.suffix\n";
	print "RMSD cutoff: Determines the distance threshold used to create clusters\n";
	print "allrmsd/lrmsd: Used to determine if pairwise prediction distances use all-atom RMSD or just ligand RMSD\n";
	print "nosizes/sizes: Option that specifies if the sizes of the different clusters should be output\n";
	print "execute/print: The former will trigger the execution of the stat generation program in the current machine, while the second ";
	print "option will only output the string representation of the command that needs to be invoked, in case the user wants to execute them ";
	print "on different machines\n";
	print "\nExample multilzerd_pairwise_cluster.pl ./predictions/ 1TB4 A,B,C,D clustersuffix 5 allrmsd sizes print\n";
}
else
{
	$directory = $ARGV[0];
	$pdbid = $ARGV[1];
	@chains = split(/,/, $ARGV[2]);
	$suffix = $ARGV[3];
	$cutoff = $ARGV[4];
	$rmsdtype = $ARGV[5];
	$printsizes = ($ARGV[6] eq "sizes");
	$mode = $ARGV[7];

	$clusterprogram = "lzerd_cluster";

	foreach $rec_chain_index (0 .. $#chains)
	{
		foreach $lig_chain_index ($rec_chain_index + 1 .. $#chains)
		{
		
			$rec_chain = $chains[$rec_chain_index];
			$lig_chain = $chains[$lig_chain_index];

			$rec_pdb = "$directory/$rec_chain-$pdbid.pdb";
			$lig_pdb = "$directory/$lig_chain-$pdbid.pdb";

			$predictionsfile = "$directory/$rec_chain-$lig_chain.out";
			
			$outputfile = "$directory/$rec_chain-$lig_chain.$suffix";

			# in case the cluster sizes are printed out, this would be the file
			$sizesfile = "$directory/$rec_chain-$lig_chain.$suffix.sizes";


			$command = "(nohup $clusterprogram -R $rec_pdb -L $lig_pdb --input $predictionsfile --cutoff $cutoff ";

			# provide the lrmsd flag only if it's requested
			if($rmsdtype eq "lrmsd")
			{
				$command .= "--lrmsd ";
			}

			# similarly, we only add the printsizes parameter if we need it
			if($printsizes)
			{
				$command .= "--printsizes $sizesfile ";
			}

			# at this point the command is complete, but we need to add the redirection
			$command .= "> $outputfile) > $outputfile.nohup 2>&1&";

			if($mode eq "execute")
			{
				system($command);
			}
			else
			{
				print "$command\n";
			}
		}
	}
	if($mode eq "execute")
	{
		print "\nProcesses are running in the background. The results will be output to files of the form <chainID>-<chainID>.<suffix> ";
		print "for each pairwise combination.\n";
	}
}
