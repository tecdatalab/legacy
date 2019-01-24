#!/usr/bin/perl

# Receives a file that contains correlations between patches
# and outputs a one column stream with the distance between geometric
# centers of patch1 and patch2.

if ($#ARGV + 1 != 1) {
	print "Usage: patch_distances.pl <correlations file>\n";
} else {
	$correlations_filename = $ARGV[0];
	open CORRFILE, $correlations_filename or die $!;

	for $corr_line(<CORRFILE>) {
		chomp($corr_line);
		@line_parts = split(/\s+/, $corr_line);
		$patch1_file = $line_parts[0];
		$patch2_file = $line_parts[1];

		$patch1_file =~ s/\.inv//;
		$patch2_file =~ s/\.inv//;

		print `../centroid_distance.py $patch1_file $patch2_file`;
	}
}
