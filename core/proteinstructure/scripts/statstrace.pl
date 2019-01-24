#!/usr/bin/perl

# When Multi-LZerD's --detailed option is used, we can use the partial results output to analyze the behavior of the best scoring elements
# as the number of generations increase
# This script will simply integrate a subset of the results held in all those files. For example, if results per generation = 1, it will generate
# a statistics file with the information corresponding to the top prediction in each of the generations

if($#ARGV + 1 != 6)
{
	print "Usage: statstrace.pl <filename suffix> <ga results directory> <template directory> <results per generation> <start generation #> <end generation #>\n";
}
else
{
	$suffix = $ARGV[0];
	$directory = $ARGV[1];
	$templates = $ARGV[2];
	$pergeneration = $ARGV[3];
	$generationindex = $ARGV[4] - 1; # to make it start from 0
	$lastgeneration = $ARGV[5];
	$headerprinted = 0;

	while($generationindex++ < $lastgeneration)
	{
		# file indices start from 1 (not from zero) so the previous increment is correct
		$filenumber = sprintf("%05d", $generationindex);
		$filename = "$directory/$filenumber-$suffix";
		$statsfile = "$filename.statstmp";
		# generate a tmp file that holds the result of calling mdga_stats
		`mdga_stats --input $filename -d $templates -n $pergeneration > $statsfile`;

		# open it and print to std output $pergeneration elements per file
		open(STATSTMP, $statsfile);
		$header = <STATSTMP>;
		
		# print the header if it's the first time
		if(!$headerprinted)
		{
			print "Generation,";
			print $header;
			$headerprinted = 1;
		}
		#print the first $pergeneration lines
		for($i = 0; $i < $pergeneration; $i++)
		{
			print "$generationindex,";
			print <STATSTMP>;
		}

		close(STATSTMP);

		#delete the tmp file
		`rm $statsfile`;
	}
}
