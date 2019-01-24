#!/usr/bin/perl

# Script: create_pdb_latest_basic.pl
# Description:
# Author: YI XIONG
# email: xiongyi128@gmail.com
# Date: Nov.11.13

if ( $#ARGV + 1 != 3 ) {
	print
"Usage: create_pdb_latest_basic.pl <inv directory> <pdb list file> <output filename>\n";
}
else {

	$inv_dir             = $ARGV[0];
	$pdb_list_file    = $ARGV[1];
	$output_filename     = $ARGV[2];
	$output_tmp_filename = "$output_filename.tmp";

	$inv_counter = 0;

	open( LIST, $pdb_list_file );
	@pdb_list_array = <LIST>;
	close(LIST);

# open the temprorary file that will contain the 3DZD in each inv without the header of the file that will be added later
	open( TMP_OUT, ">$output_tmp_filename" );
	foreach $pdb_id (@pdb_list_array) {

		chomp($pdb_id);
		$inv_filename = $pdb_id . ".pdb.inv";
		$file         = "$inv_dir/$inv_filename";

		if ( -e $file ) {
			
			print TMP_OUT $inv_filename;
			open( IN, $file );
			# read the information in the file and write it
			while (<IN> ) {
				chomp($_);
				print TMP_OUT " $_";
			}
			close(IN);
			print TMP_OUT "\n";
			$inv_counter++;
		}
	}

	# close tmp files
	close(TMP_OUT);

	# create the final file and write the header
	open( OUT, ">$output_filename" );
	print OUT "ZERNIKE $inv_counter\n\n";
	close(OUT);

	# concatenate the tmp file to the final one
	`cat $output_tmp_filename >> $output_filename`;

	#remove the temporary files
	`rm $output_tmp_filename`;
}

print "Success\n";
