#!/usr/bin/perl

# INV files contain the zernike descriptors that describe the shape of proteins
# This program reads the two descriptor files and calculates the euclidean distance between
# both sets of values

if($#ARGV + 1 != 2)
{
	print "Usage: invdistance.pl <inv file> <inv file>\n";
}
else
{
	$inv1 = $ARGV[0];
	$inv2 = $ARGV[1];

	#open both files
	open(INV1, $inv1);
	open(INV2, $inv2);
	# the first line contains the number of descriptors
	$inv1_total = <INV1>;
	$inv2_total = <INV2>;
	chomp($inv1_total);
	chomp($inv2_total);
	if($inv1_total != $inv2_total) # this is an error, they should be the same length
	{
		print "The inv files have a different number of Zernike Descriptors\n";
	}
	else
	{
		$squared_distance = 0;
		for($index = 0; $index < $inv1_total; $index++)
		{
			$descriptor1 = <INV1>;
			chomp($descriptor1);
			$descriptor2 = <INV2>;
			chomp($descriptor2);

			$squared_distance += ($descriptor1 - $descriptor2) ** 2;

		}
		$distance = sqrt($squared_distance);

		# print the distance value
		print "$distance\n";
	}
	close(INV1);
	close(INV2);
}
