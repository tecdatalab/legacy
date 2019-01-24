#!/usr/bin/perl

use strict;
use warnings;

sub download_map {
  my ($local_file, $ftp_file) = @_;
  our $wget_logfile;
  my $wget_result = system("wget -a $wget_logfile -O $local_file -- ftp.ebi.ac.uk$ftp_file");
  if ($wget_result == 0) {
    `gunzip $local_file`;
  } else {
    print "Error downloading map $ftp_file\n";
    `rm -f $local_file`;
  }
}

sub download_xml {
  my ($local_file, $ftp_file) = @_;
  our $wget_logfile;
  my $wget_result = system("wget -a $wget_logfile -O $local_file -- ftp.ebi.ac.uk$ftp_file");
  if ($wget_result != 0) {
    print "Error downloading XML header $ftp_file\n";
    `rm -f $local_file`;
  }
}

sub download_image {
  my ($tmp_image_file, $local_file, $ftp_file) = @_;
  our $wget_logfile;
  my $wget_result = system("wget -a $wget_logfile -O $tmp_image_file -- ftp.ebi.ac.uk$ftp_file");
  if ($wget_result == 0) {
    `convert $tmp_image_file $local_file`;
  } else {
    print "Error downloading image $ftp_file\n";
  }
  `rm -f $tmp_image_file`;
}

# Main entry point
if ($#ARGV != 2) {
  die "Usage: emdb_download_updated_files.pl <old timestamps> <new timestamps> <wget logfile>\n";
}

use Date::Parse;

my %old_timestamps = ();
my ($old_file, $new_file);
open($old_file, $ARGV[0]) || die "Could not open ". $ARGV[0] . "\n";

my $is_header_processed = 0;

for my $old_line(<$old_file>) {
  if (!$is_header_processed) {
    $is_header_processed = 1;
    next;
  }
  chomp($old_line);
  my ($id, $mapdate, $xmldate, $imagedate, $mapfile, $xmlfile, $imagefile) =
      split(/\t/, $old_line);
  push(@{$old_timestamps{$id}}, (str2time($mapdate), str2time($xmldate),
      str2time($imagedate), $mapfile, $xmlfile, $imagefile));
}

close($old_file);

open($new_file, $ARGV[1]) || die "Could not open ". $ARGV[1] . "\n";
$is_header_processed = 0;

our $wget_logfile = $ARGV[2];
`rm -f $wget_logfile`;

# First clean the directories used to download files
`rm -Rf img xml map`;
`mkdir img`;
`mkdir xml`;
`mkdir map`;

for my $new_line(<$new_file>) {
  if (!$is_header_processed) {
    $is_header_processed = 1;
    next;
  }
  chomp($new_line);
  my ($id, $mapdate, $xmldate, $imagedate, $mapfile, $xmlfile, $imagefile) =
      split(/\t/, $new_line);
  my $downloadedmap = "map/" . $id . ".map.gz";
  my $downloadedxml = "xml/emd-" . $id . ".xml";
  my $downloadedimage = "img/emd-" . substr($imagefile, rindex($imagefile, "/") + 1);

  if (not exists $old_timestamps{$id}) {
    # Since it didn't exist, we need to update it.
    download_map($downloadedmap, $mapfile);
    download_xml($downloadedxml, $xmlfile);
    download_image($downloadedimage, "img/$id.gif", $imagefile);
  } else {
    my @old_dates = @{$old_timestamps{$id}};
    if ($old_dates[0] < str2time($mapdate) ||
    	$old_dates[1] < str2time($xmldate) ||
    	$old_dates[2] < str2time($imagedate)) {
      download_map($downloadedmap, $mapfile);
      download_xml($downloadedxml, $xmlfile);
      download_image($downloadedimage, "img/$id.gif", $imagefile);
    }
  }
}
close($new_file);
