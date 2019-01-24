#!/usr/bin/perl

my $num           = $ARGV[0];
my $volume_filter = $ARGV[1];
my $query_id_list = $ARGV[2];
my $inv_dir       = $ARGV[3];
my $query_dir       = $ARGV[4];
my $db_prog       = "querybatchsearchidlist.pl";

$query_id_list =~ s/\A\s+|\s+\z//g;
@listarray = split( /,/, $query_id_list );

`rm $query_dir.zip`;
`rm -r $query_dir`;
`mkdir $query_dir`;

if ( $volume_filter eq "-nofilter" ) {

	foreach my $i ( 0 .. $#listarray ) {

		$query_id    = $listarray[$i];
		$queryresult = "EMD-" . $query_id . ".hit";

		open( OUT, ">$query_dir/$queryresult" );

		print OUT "**********************************\n";
		print OUT "Results for query: EMD-$query_id\n";
		print OUT "**********************************\n";
		print OUT "RANK\tEMDB_ID\t\tEUC_DIST\n";
		print OUT "**********************************\n";

		for (`./$db_prog $query_id $volume_filter $inv_dir $num`) {
			next if !$_;
			chomp($_);
			( $rank, $template_id, $dist ) = split( /\t/, $_ );
			print OUT $rank 
			  . "\tEMD-"
			  . $template_id . "\t"
			  . $dist . "\n";
		}

		print OUT "**********************************\n";
		close(OUT);
	}
}
elsif ( $volume_filter eq "-filter" ) {
	foreach my $i ( 0 .. $#listarray ) {

		$query_id    = $listarray[$i];
		$queryresult = "EMD-" . $query_id . ".hit";

		open( OUT, ">$query_dir/$queryresult" );

		print OUT "*************************************************\n";
		print OUT "Results for query: EMD-$query_id\n";
		print OUT "*************************************************\n";
		print OUT "RANK\tEMDB_ID\t\tVOLUME_RATIO\tEUC_DIST\n";
		print OUT "*************************************************\n";

		for (`./$db_prog $query_id $volume_filter $inv_dir $num`) {
			next if !$_;
			chomp($_);
			( $rank, $template_id, $dist, $ratio ) = split( /\t/, $_ );
			print OUT $rank 
			  . "\tEMD-"
			  . $template_id . "\t"
			  . $ratio . "\t\t"
			  . $dist . "\n";
		}

		print OUT "*************************************************\n";
		close(OUT);
	}
}

