#!/usr/bin/perl

#
# Script: perl query.pl id_query
#
# Description:
# Author: Yi Xiong
# email: xiongyi128@gmail.com
# Date: Dec.11.2013

sub compute_distance {
	my ( $file1, $file2 ) = @_;
	open FILE1, $file1 or die $!;
	open FILE2, $file2 or die $!;
	my @lines1 = <FILE1>;
	my @lines2 = <FILE2>;
	close FILE2;
	close FILE1;
	chomp(@lines1);
	chomp(@lines2);

	my $sum = 0;

	for (
		my $descriptor_index = 1 ;
		$descriptor_index <= $#lines1 ;
		$descriptor_index++
	  )
	{
		$sum += ( $lines1[$descriptor_index] - $lines2[$descriptor_index] )**2;
	}

	my $dist = sqrt($sum);

	$dist;
}

sub is_filtered_by_resolution {
	my($current_resolutions_hash, $emdb_id, $min_param, $max_param) = @_;
	my $query_resolution = -1;
	if (exists $current_resolutions_hash->{$emdb_id}) {
		$query_resolution = $current_resolutions_hash->{$emdb_id};
	}
	if ($min_param == -1 && $max_param == -1) {
		return 0;
	} elsif ($query_resolution == -1) {
		# i.e. if at least one of the filters was provided but the current EM map
		# doesn't have a resolution assigned, we want to ignore it to just return
		# maps that strictly comply with the filters
		return 1;
	} elsif (($min_param > 0 && $query_resolution < $min_param) ||
		 ($max_param > 0 && $query_resolution > $max_param)) {
		 return 1;
	} else {
		return 0;
	}
}

$id_query      = $ARGV[0];
$volume_filter = $ARGV[1];
$inv_dir = $ARGV[2];
$min_resolution = $ARGV[3];
$max_resolution = $ARGV[4];
if ($#ARGV == 6 && $ARGV[5] eq "-upload") {
	$search_upload = 1;
	$query_volume = $ARGV[6];
} else {
	$search_upload = 0;
}

# Read volumes
my %volumehash;
open( IN, "../db/volume.txt" );
@inarray = <IN>;
close(IN);
chomp(@inarray);

foreach $line (@inarray) {
	my @tmparray = split( /\s+/, $line );
	if ( $tmparray[1] > 0 ) {
		$volumehash{ $tmparray[0] } = $tmparray[1];
	}
}

# Read resolutions
my %resolutionhash;
open( INRESOLUTIONS, "../db/resolutions.txt" );
@resolutionsarray = <INRESOLUTIONS>;
close(INRESOLUTIONS);
chomp(@resolutionsarray);

foreach $line (@resolutionsarray) {
	my @tmparray = split( /\s+/, $line );
	if ( $tmparray[1] > 0 ) {
		$resolutionhash{ $tmparray[0] } = $tmparray[1];
	}
}

opendir( DIR, $inv_dir );
@filearray = grep( /\.inv$/, readdir(DIR) );
closedir(DIR);

my %iddisthash;
foreach $filename (@filearray) {
	$id_template = substr( $filename, 0, 4 );
	if ( $id_query ne $id_template || $search_upload) {
		my $filename1 = "$inv_dir/$id_query.inv";
		my $filename2 = "$inv_dir/$id_template.inv";
		if ($search_upload) { # Overwrite the default path
			$filename1 = $id_query;
		}
		my $dist = compute_distance( $filename1, $filename2 );
		$dist = sprintf "%.3f", $dist;
		if ( $dist >= 0 ) {
			$iddisthash{$id_template} = $dist;
		}
	}
}

if ( $volume_filter eq "-nofilter" ) {
	my $num = 0;
	foreach
	  my $key ( sort { $iddisthash{$a} <=> $iddisthash{$b} } keys %iddisthash )
	{
		if (!is_filtered_by_resolution(\%resolutionhash, $key, $min_resolution, $max_resolution)) {
			++$num;
			if ( $num <= 20 ) {
				print "$key\t$iddisthash{$key}\t$resolutionhash{$key}\n";
			}
		}
	}
}
elsif ( $volume_filter eq "-filter" ) {
	my $num = 0;
	if (!$search_upload) {
		$query_volume = $volumehash{$id_query};
	} # else it should have been passed as parameter
	foreach
	  my $key ( sort { $iddisthash{$a} <=> $iddisthash{$b} } keys %iddisthash )
	{
		if ( exists $volumehash{$key} ) {
			$ratio = $volumehash{$key} / $query_volume;
			if ( $ratio >= 0.8 && $ratio <= 1.2 ) {
				if (!is_filtered_by_resolution(\%resolutionhash, $key,
							$min_resolution, $max_resolution)) {
					++$num;
					if ( $num <= 20 ) {
						print "$key\t$iddisthash{$key}\t$resolutionhash{$key}\n";
					}
				}
			}
		}
	}
}
