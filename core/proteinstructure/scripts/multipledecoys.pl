#!/usr/bin/perl

# Before analyzing multiple docking results it is desirable to know if the pair wise
# predictions generated by LZerD are good or not, to create a complete multi-chain predicted complex
# To test this we can create sample .ga.out files that define first, the connections between chains and
# second, the prediction numbers for each connection (LZerD original rank)
# In order to know which pair wise predictions are good we can use the CSTATS_VDOCK program and select
# the predictions that have good RMSD and/or fnat values
# With this information we can manually create "summary files" (processed by this script) which
# contain a subset of the standard .ga.out header first, and then they specify the rank of several
# predictions for each pair of chains. We only need to specify a subset of all connections, specifically
# a subset of the connections that creates a spanning tree to connect all chains
#
# Finally, this program creates a .ga.out file that combines all the different predictions numbers for each connection
# and then we can feed that file to the multiple docking statistics program

# Recursive function that helps in the recursive creation of prediction sets
sub printpredictionset
{
        my ($prefix, $refpairwisepredictions, $pairwiseindex) = @_;
	my @pairwisepredictions = @$refpairwisepredictions;
	if($pairwiseindex > $#pairwisepredictions)
	{
		# just print what we have so far with dummy score values
		push(@finalresult, "$prefix\t0\t0\t0\n");
	}
	else
	{
		# get the next pairwise set and for each prediction recursively call the function
		my $predref = $pairwisepredictions[$pairwiseindex];
		my @pred = @$predref;
		my $predindex = 0;
		while($predindex < $#pred + 1)
		{
			my $value = $prefix . $pred[$predindex];
			if($pairwiseindex < $#pairwisepredictions)
			{
				$value .= ";";
			}
			printpredictionset($value, \@pairwisepredictions, $pairwiseindex + 1);
			$predindex++;
		}
	}
}

if ($#ARGV != 0)
{
	print "Usage: multipledecoys.pl <summary file>\n";
}
else
{
	@finalresult = ();
	open(SUMMARY, "$ARGV[0]");
	@lines = <SUMMARY>;
	chomp(@lines);
	
	# it is expected that the file starts with a header similar to the one found in GA output files
	# headernames are ignored (it is expected that the order is the following
	($headername, $pdbid) = split("\t", $lines[0]);
	($headername, $allchains) = split("\t", $lines[1]);
	($headername, $clashes) = split("\t", $lines[2]);
	($headername, $scoretype) = split("\t", $lines[3]);
	($headername, $decoyprogram) = split("\t", $lines[4]);

	# associate each of the chain id's to the position they were found in
	%chains;
	@chainids = split(",", $allchains);
	for($i = 0; $i < $#chainids + 1; $i++)
	{
		$chains{$chainids[$i]} = $i;
	}


	# the fifth line should be empty and starting from line 6 we should have <chainid>-<chainid>
	# followed by several lines with 4 numbers: the pairwise prediction number, iRMSD, LRMSD and fnat
	# to create a sample results file using those predictions we just need the prediction number
	# (the other measures are there for clarity purposes, if someone reads the file)
	$lineindex = 6;
	$pairwiseindex = 0;
	@pairwisepredictions;
	while($lineindex < $#lines + 1)
	{
		my @onepairwiseset = ();
		# the first line of a segment should be the chainid-chainid pair
		@pairofchains = split("-", $lines[$lineindex++]);

		# extract the chain id index, acoording to the order in which they were provided
		$receptorindex = $chains{$pairofchains[0]};
		$ligandindex = $chains{$pairofchains[1]};

		# at this point it's over the first prediction
		# extract predictions until the lines finish or until another pairwise set comes up 
		while(($lineindex < $#lines + 1) && ($lines[$lineindex] ne ""))
		{
			($prednumber,$irmsd,$lrmsd,$fnat) = split(/\s+/, $lines[$lineindex++]);
			push(@onepairwiseset, "$receptorindex,$ligandindex,$prednumber");
		}
		# when it finishes the current set of pairwise predictions, we add them to the complete set
		$reftoarray = \@onepairwiseset;
		push(@pairwisepredictions, \@onepairwiseset);
		
		$reftoarray = $pairwisepredictions[0];
		# this jumps over the empty line
		$lineindex++;

	}


	printpredictionset("", \@pairwisepredictions, 0);
	
	$popsize = $#finalresult + 1;
	print "PDBID\t$pdbid\nChains\t$allchains\nGenerations\t1\nPopulationSize\t$popsize\n";
	print "ClashThreshold\t$clashes\nScoreType\t$scoretype\nDecoyProgram\t$decoyprogram\n\n";

	print @finalresult;

	close(SUMMARY);
}
