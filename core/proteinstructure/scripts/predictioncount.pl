#!/usr/bin/perl

# If Multi-LZerD is invoked with the --detailed option it will generate an output file for each generation
# This program analyzes all of those output files in order to create a report of the number of times
# one particular edge is used, for example, chain A joined to chain B using prediction #456
#
# Since all those files have the same base name, just changing the generation number, we can generate results
# by providing the suffix of the filename, the directory where they are stored, and the range of generations
# (i.e. from generation 5 to 250) that we want to analyze.
# The final parameter, step size, when assigned a value of 1 will instruct the script to generate a count report
# separately for each generation but, since we can have a large number of generations, we may want to generate
# count reports every several generations. For example, a 5 would instruct the script to integrate the counts for 
# every set of 5 generations

sub getrangecount
{
	my ($suffix, $directory, $generationindex, $lastgeneration) = @_;

	# this is the hash table that stores the number of times that particular edges occur (initially empty)
	my %edgecount = ();

	while($generationindex++ < $lastgeneration)
	{
		# file indices start from 1 (not from zero) so the previous increment is correct
		my $filenumber = sprintf("%05d", $generationindex);
		my $filename = "$directory/$filenumber-$suffix";

		# open the output file corresponding to the current generation
		open(GAFILE, $filename);
		my $pdbidheader = <GAFILE>;
		my $chainsheader = <GAFILE>;
		chomp($chainsheader);

		# extract the chain names so we can use them to "translate" the chain numbers
		@chains = split(/\t+|s+|,/, $chainsheader);
		# this removes the header name and the tab and leaves just the chains
		shift(@chains);
		shift(@chains);

		# read the rest of the header lines but ignore them
		my $headersfinished = 0;
		while(!$headersfinished)
		{
			# an empty line marks the beginning of the actual data
			my $ignoredheader = <GAFILE>;
			chomp($ignoredheader);
			if($ignoredheader eq '')
			{
				$headersfinished = 1;
			}
		}

		# now we read each prediction and increment the appropriate counters
		while(<GAFILE>) {
			my $nextline = <GAFILE>;
			chomp($nextline);
			# each prediction line is composed of the actual edge description, pairwise score, clashes and physics score
			my @predictionparts = split(/\t+/, $nextline);
			# we only use the edge description and separate each edge description (separated by ;)
			my @edges = split(/;/,$predictionparts[0]);

			#go through each of the edges and extract the receptor, ligand and prediction number (comma separated)
			for my $edge(@edges)
			{
				my @edgeparts = split(/,/, $edge);
				my $receptor = $edgeparts[0];
				my $ligand = $edgeparts[1];
				my $predictionnumber = $edgeparts[2];

				$edgecount{$receptor}{$ligand}{$predictionnumber}++;
			}
		}

		close(GAFILE);
	}

	# return the edge count
	return %edgecount;

}

if($#ARGV + 1 != 5)
{
	print "Usage: predictioncount.pl <filename suffix> <ga results directory> <start generation #> <end generation #> <step size>\n";
}
else
{
	$suffixparam = $ARGV[0];
	$directoryparam = $ARGV[1];
	$generationindexparam = $ARGV[2] - 1; # to make it start from 0
	$lastgenerationparam = $ARGV[3];
	$stepsizeparam = $ARGV[4];

	# go from "generation index" to "last generation" analyzing "step size" generations at a time
	
	# keep a start and end position for each range
	$start = $generationindexparam;
	$end = $start + $stepsizeparam;

	# create a tmp file that will store the counts prior to sorting
	$presortfilename = "tmp$suffixparam.nosort";

	# the following array will store all the tags used to identify the ranges
	@tags = ();

	open(RESULTSFILE, ">$presortfilename");

	while($start < $lastgenerationparam)
	{
		# count the elements for the current range and store the hash table with the results
		%currentedgecount = getrangecount($suffixparam, $directoryparam, $start, $end);
		
		#These are used to identify the generation range that we are dealing with at the moment
		$starttag = sprintf("%05d", $start + 1);
		$endtag = sprintf("%05d", $end);
		$tag = "$starttag-$endtag";
		push(@tags, $tag);

		# output the results corresponding to the current range
		for $receptorkey (keys(%currentedgecount))
		{
			for $ligandkey (keys(%{$currentedgecount{$receptorkey}}))
			{
				for $predictionkey (keys(%{$currentedgecount{$receptorkey}{$ligandkey}}))
				{
					$count = $currentedgecount{$receptorkey}{$ligandkey}{$predictionkey};
					$receptorname = $chains[$receptorkey];
					$ligandname = $chains[$ligandkey];
					print RESULTSFILE "$receptorname,$ligandname,$predictionkey,$tag,$count\n";
				}
			}
		}


		#update start and end
		$start += $stepsizeparam;
		$end = $start + $stepsizeparam;
	}

	close(RESULTSFILE);

	
	# we have a results file up to this point. We want to sort it for convenience, first, and second
	# because we also want to make sure that for each edge that we found, there are entries that describe
	# the number of times it happens in each range.
	# For example, edge A-B,400 could have happened just in the second range and then it didn't happen again,
	# however, we also want to make sure that we have a row in the output that says that in the first range
	# it didn't occur (namely, add a zero)
	
	# first print an initial header line
	print "Receptor,Ligand,Prediction";
	for $tagtitle(@tags)
	{
		print ",$tagtitle";
	}
	
	# if we have N tags, we expect that a sequence of rows will have tags 0, 1, 2...N
	# and then repeat 0, 1, 2....N (we will cycle through them and fill the gaps
	$expectedtagindex = 0;
	@sortedentries = `sort -t , -k 1,2 -k 3,3n -k 4,4 $presortfilename`;

	#use these variables to control if at each iteration we are analyzing the same edge, but in a different range
	$previousreceptor = "";
	$previousligand = "";
	$previousprediction = -1;


	for($entryindex = 0; $entryindex <= $#sortedentries;)
	{
		$sortedentry = $sortedentries[$entryindex];
		($receptorname,$ligandname,$predictionkey,$tag,$count) = $sortedentry =~ /(.*),(.*),(.*),(.*),(.*)/;

		# first check if we are beginning a new row either because the expected tag is zero or because
		# the receptor/ligand/prediction are different from the previous
		
		$isdifferent = $previousreceptor ne $receptorname
					|| $previousligand ne $ligandname || $previousprediction ne $predictionkey;
		if($expectedtagindex == 0)
		{
			# print the "header section" (i.e. the initial columns that are not counts)
			print "\n$receptorname,$ligandname,$predictionkey";
		}

		# and now print the count corresponding to tags
		$expectedtag = $tags[$expectedtagindex];

		# if we are the beginning and the tag matches we just print the count
		# also, if the tag matches and we are still on the same "prediction" we print the value
		if($expectedtag eq $tag && ($expectedtagindex == 0 || !$isdifferent))
		{
			# then just print the count and update the expected tag index
			print ",$count";
			$expectedtagindex = $expectedtagindex < $#tags ? $expectedtagindex + 1 : 0;

			#update the entry index to the next (in the following cases this doesn't always happen)
			$entryindex++;
			# update all the "previous" variables
			$previousreceptor = $receptorname;
			$previousligand = $ligandname;
			$previousprediction = $predictionkey;
		}
		else
		{
			# there is some type of gap
			# if the expected tag is zero at this point, it means that the tag wasn't equal
			if($expectedtagindex == 0)
			{
				if($isdifferent)
				{
					# this means that we need to fill on the left part
					while($tag ne $expectedtag)
					{
						$expectedtag = $tags[++$expectedtagindex];
						print ",0";
					}
					# update all the "previous" variables
					$previousreceptor = $receptorname;
					$previousligand = $ligandname;
					$previousprediction = $predictionkey;
				}
				else # this shouldn't happen, abort execution
				{
					die "Error: Tag expected was zero but program is still on the same entry\n";
				}
			}
			elsif($isdifferent) # but note that the expected tag is not zero
			{
				# we need to right-fill until we run out of tags
				while($expectedtagindex <= $#tags)
				{
					$expectedtagindex++;
					print ",0";
				}
				# at the end reset the expectet to 0
				$expectedtagindex = 0;
			}
			elsif($expectedtag ne $tag)
			{
				# we need to fill a gap in the middle
				while($tag ne $expectedtag)
				{
					$expectedtag = $tags[++$expectedtagindex];
					print ",0";
				}
			}
			else
			{
				print "\n";
				print $expectedtagindex == 0;
				print "\n";
				print $isdifferent;
				print "\n";
				die "Error: Unexpected failure\n";
			}
		}
	}

	#make sure that the last line is complete
	if($expectedtagindex != 0)
	{
		while($expectedtagindex <= $#tags)
		{
			$expectedtagindex++;
			print ",0";
		}
	}

	print "\n";

	`rm $presortfilename`;
}
