#!/usr/bin/perl

# Script: perl upload_to_inv.pl file.map contour representation output.inv

use File::Temp qw/ tempfile tempdir /;

sub compute_one_inv_set {
	my($source, $current_contour) = @_;
	#my($tmp_handle, $tmp_filename_full) = tempfile(UNLINK => 1, SUFFIX => '.inv');
	#my($tmp_handle, $tmp_filename) = tempfile(UNLINK => 1);
	my $tmp_filename = $source . "_tmp";
	`cp $source $tmp_filename.map`;
	`uploadscripts/map2zernike $tmp_filename.map -n 20 -c $current_contour -p $tmp_filename`;
	#`uploadscripts/em_to_3dzd.sh $source $current_contour $tmp_filename`;
	my ($length) = `wc $tmp_filename.inv` =~ m/(\d+)/;
	$length--; # Because it includes a one-line header with the size
	my $descriptors = `tail -n $length $tmp_filename.inv`;
	`rm $tmp_filename.inv $tmp_filename.3dzm $tmp_filename.map`;
	($length, $descriptors);
}

if ($#ARGV != 3) {
	print("Usage: upload_to_inv.pl file.map contour representation output.inv\n");
} else {
	$input = $ARGV[0];
	$contour = $ARGV[1];
	$representation = $ARGV[2];
	$output = $ARGV[3];

	# First get the min and max densities to determine the 1/3 and 2/3 thresholds
	$volume_result = `uploadscripts/em_volume $contour $input`;
	($min_density, $max_density) = $volume_result =~ m/Min\/Max densities: (.*) (.*)\s/;
	
	$density_range = ($max_density - $contour);
	$one_third_contour = $contour + (1/3) * $density_range;
	$two_thirds_contour = $contour + (2/3) * $density_range;

	# Recommended contour is computed always
	($total_descriptors, $all_descriptors) = compute_one_inv_set($input, $contour);

	if ($representation eq "recommendonethird" ||
		$representation eq "recommendonethirdtwothird") {
		($current_length, $current_descriptors) =
			compute_one_inv_set($input, $one_third_contour);
		$total_descriptors += $current_length;
		$all_descriptors .= $current_descriptors;
	}

	if ($representation eq "recommendtwothird" ||
		$representation eq "recommendonethirdtwothird") {
		($current_length, $current_descriptors) =
			compute_one_inv_set($input, $two_thirds_contour);
		$total_descriptors += $current_length;
		$all_descriptors .= $current_descriptors;
	}
	
	# Since we don't have the XML with the std dev then assume stddev of 1.
	if ($representation eq "recommendonestd") {
		($current_length, $current_descriptors) =
			compute_one_inv_set($input, $contour + 1);
		$total_descriptors += $current_length;
		$all_descriptors .= $current_descriptors;
	}
	# Write to the output file
	open(OUTPUT_HANDLE, ">", $output);

	print OUTPUT_HANDLE "$total_descriptors\n$all_descriptors";

	close(OUTPUT_HANDLE);
}
