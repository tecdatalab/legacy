#!/usr/bin/perl
if($#ARGV + 1 != 1)
{
	print "Usage: score_lzerd.pl <lzerd output>\n";
}
else
{
	$lzerdfile = @ARGV[0];
	# Weights used when considering all 12 terms
	@weights12 = (53.074, 151.294, 0.15, 1.373, 13.235, -0.097, 15.93, -0.734, 0, 0, 27.635, 79.063);
	# Weights used when using 9 terms, excluding overall van der waals, electrostatics and solvation
	@weights9 = (0, 81.9602, 0.1916921, 0, 10.57919, -0.05739181, 12.37243, -0.5037819, 0, 2.853197, 0, 58.7541);

	open(LZERD_FILE, $lzerdfile);
	open(OUTPUT12, ">$lzerdfile.tmp.12");
	open(OUTPUT9, ">$lzerdfile.tmp.9");
	
	$index = 1;
	while(<LZERD_FILE>)
	{
		chomp($_);
		@terms = split(/ +/, $_);
		$score12 = 0;
		$score9 = 0;
		for($i=0; $i < 12; $i++)
		{
			$score12 = $score12 + ($weights12[$i] * $terms[$i]);
			$score9 = $score9 + ($weights9[$i] * $terms[$i]);
		}
		
		print OUTPUT12 "$index\t$score12\n";
		print OUTPUT9 "$index\t$score9\n";

		$index++;
	}


	close(OUTPUT12);
	close(OUTPUT9);
	close(LZERD_FILE);

	# and now sort them to create the two rankings
	`sort -n -k 2,2 $lzerdfile.tmp.12 > $lzerdfile.rank12`;
	`sort -n -k 2,2 $lzerdfile.tmp.9 > $lzerdfile.rank9`;

	#remove temp files
	`rm -f $lzerdfile.tmp.12 $lzerdfile.tmp.9`;

}
