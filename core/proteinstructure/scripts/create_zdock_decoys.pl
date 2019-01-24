#!/usr/bin/perl
if($#ARGV + 1 != 2)
{
	print "This command assumes that there is a separate file for each chain and that its name follows the <chainID>-<PDB ID>.pdb.ms format\n";
	print "The script will begin the execution of zdock for each pairwise combination and leave processes running in the background.\n";
	print "\nUsage: create_zdock_decoys <PDB ID> <comma separated chain names>\n";
	print "Example create_zdock_decoys 1Z5S A,B,C,D\n";
}
else
{
	$chain_id = $ARGV[0];
	@chains = split(/,/, $ARGV[1]);
	foreach $rec_chain_index (0 .. $#chains)
	{
		foreach $lig_chain_index ($rec_chain_index + 1 .. $#chains)
		{
		
			$rec_chain = $chains[$rec_chain_index];
			$lig_chain = $chains[$lig_chain_index];
			$rec_file = "$rec_chain-$chain_id.pdb.ms";
			$lig_file = "$lig_chain-$chain_id.pdb.ms";
			$out_file = "$rec_chain-$lig_chain.out";
			$nohup_file = "$rec_chain-$lig_chain.nohup";

			system("nohup ~/bin/zdock -N 54000 -D -o $out_file -R $rec_file -L $lig_file &> $nohup_file &");
		}
	}
	print "Processes are running in the background. The results will be output to files of the form <chainID>-<chainID>.out ";
	print "for each pairwise combination.\nA series of .nohup files will be generated too; errors generated by zdock can be seen there.\n";
}