#!/usr/bin/perl

# In order to test how good pairwise LZerD predictions are, RMSD and fnat stats can be generated.
# This script assumes we're using bound docking and that pairwise predictions have already been generated
# and reduced to a manageable number. It invokes a program that will calculate the LRMSD, iRMSD and fnat
# between the predictions and the native conformation

if($#ARGV + 1 != 5)
{
	print "Assuming pairwise LZerD predictions have been generated already, and that they follow the <chain>-<chain>.out naming format, ";
	print "this script will generate the LRMSD, iRMSD and fnat statistics for each of the predictions provided in the LZerD output file.\n";
	print "\nUsage: multilzerd_pairwise_stats <directory> <PDB ID> <comma separated chain names> <output suffix> <execute/print>\n\n";
	print "directory: Where the LZerD output files are located\n";
	print "PDB ID and chain names: Used to infer the names of PDB and LZerD output files in the directory\n";
	print "Output suffix: Files generated will have names like <chain>-<chain>.<suffix>\n";
	print "execute/print: The former will trigger the execution of the stat generation program in the current machine, wile the second ";
	print "option will only output the string representation of the command that needs to be invoked, in case the user wants to execute them ";
	print "on different machines\n";
	print "\nExample multilzerd_pairwise_stats.pl ./predictions/ 1TB4 A,B,C,D lzerdstats print\n";
}
else
{
	$directory = $ARGV[0];
	$pdbid = $ARGV[1];
	@chains = split(/,/, $ARGV[2]);
	$suffix = $ARGV[3];
	$mode = $ARGV[4];

	$statsprogram = "CSTATS_VDOCK";
	$lineaddprogram = "addline.sh";

	# generate a stats file for each pairwise prediction
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

			$command = "(nohup $statsprogram $rec_pdb $lig_pdb $rec_pdb $lig_pdb $predictionsfile > $outputfile 2>&1; " .
					"addline.sh $outputfile) > $outputfile.nohup 2>&1&";

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
