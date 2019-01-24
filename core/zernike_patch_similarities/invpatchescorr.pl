#!/usr/bin/perl

# Before running this script:
# 1. Run patches.py receptor.pdb.ms 20
# 2. Run patches.py ligand.pdb.ms 20
# 3. Run zernike.pl on all patches generated
# After these are done, provide the two prefixes used to name the patches
# and this script will compute the correlation coefficient between all pairs
# of patches.

sub compute_correlation {
	my ($filename1, $filename2) = @_;
	open FILE1, $filename1 or die $!;
	open FILE2, $filename2 or die $!;
	my @lines1 = <FILE1>;
	my @lines2 = <FILE2>;
	close FILE2;
	close FILE1;
	chomp(@lines1);
	chomp(@lines2);

	my $sum1 = 0;
	my $sum2 = 0;

	for(my $descriptor_index = 1; $descriptor_index <= $#lines1;
			$descriptor_index++) {
		$sum1 += $lines1[$descriptor_index];
		$sum2 += $lines2[$descriptor_index];
	}

	my $average1 = $sum1 / ($#lines1);
	my $average2 = $sum2 / ($#lines2);

	my $numerator = 0;
	my $denominator1 = 0;
	my $denominator2 = 0;

	for(my $descriptor_index = 1; $descriptor_index <= $#lines1;
			$descriptor_index++) {
		my $diff1 = ($lines1[$descriptor_index] - $average1);
		my $diff2 = ($lines2[$descriptor_index] - $average2);

		$numerator += $diff1 * $diff2;
		$denominator1 += $diff1 * $diff1;
		$denominator2 += $diff2 * $diff2;
	}
	
	my $r = $numerator / (sqrt($denominator1) * sqrt($denominator2));


	print "$filename1 $filename2 $r\n";
}

if ($#ARGV + 1 != 2) {
	print "Usage: invpatchescorr.pl <prefix name1> <prefix name2>\n";
} else {
	$prefix1 = $ARGV[0];
	$prefix2 = $ARGV[1];
	for $file1 (`ls $prefix1*atom*.pdb.ms.inv`) {
		chomp($file1);
		for $file2 (`ls $prefix2*atom*.pdb.ms.inv`) {
			chomp($file2);
			compute_correlation($file1, $file2);
		}
	}
}
