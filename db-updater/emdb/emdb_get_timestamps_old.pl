#!/usr/bin/perl

use strict;
use warnings;

use Date::Parse;
use Date::Format;
use IO::Handle;
use Net::FTP;

sub get_date_from_dir_entry {
  my ($dir_line) = @_;
  my @dir_parts = split(/\s+/, $dir_line);
  my $original_date_string = "$dir_parts[5] $dir_parts[6] $dir_parts[7]";
  my $formatted_date = time2str("%D", str2time($original_date_string));
  return $formatted_date;
}

sub get_filename_from_dir_entry {
  my ($dir_line) = @_;
  my $dir_filename = substr($dir_line, rindex($dir_line, " ") + 1);
  return $dir_filename;
}

if ($#ARGV != 0) {
  die "Usage: emdb_get_timestamps.pl <output_file>\n";
}

my $top_directory = "/pub/databases/emdb/structures/";
my $output_file;
open($output_file, ">" . $ARGV[0]) || die "Could not open ". $ARGV[0] . "\n";

STDERR->autoflush(1);
STDOUT->autoflush(1);
$output_file->autoflush(1);

print localtime() . "\n";
print "Connecting to EBI's EMDB...\n";
my $ftp = Net::FTP->new("ftp.ebi.ac.uk", Debug => 0, Passive => 1);
$ftp->login("anonymous", "") || die $ftp->message;

print "Connection established, getting timestamps...\n";

$ftp->cwd($top_directory) || die $ftp->message;

print $output_file "EMDBID\tMapModification\tXmlModification\tImageModification"
                   . "\tMapFile\tXmlFile\tImageFile\n";

my @ls_results = $ftp->ls();

$ftp->quit();

foreach my $directory (@ls_results) {
  my $emdbid = $directory;
  $emdbid =~ s/EMD-//;
  # Create a new local connection per EMDB ID. It's less efficient timewise,
  # but dropped connections have been seen frequently.

  print "Checking EMDB ID " . $emdbid . "...\n";
  $ftp = Net::FTP->new("ftp.ebi.ac.uk", Debug => 0, Passive => 1);
  $ftp->login("anonymous", "") || die $ftp->message;

  my $map_cwd = $top_directory . $directory . "/map/";
  my $xml_cwd = $top_directory . $directory . "/header/";
  my $image_cwd = $top_directory . $directory . "/images/";

  my @xml_dir = ();
  if ($ftp->cwd($xml_cwd)) {
    @xml_dir = $ftp->dir();
  } else {
    print "Error in header directory " . $emdbid . ". Message: " . $ftp->message;
  }
  my @image_dir = ();
  if ($ftp->cwd($image_cwd)) {
    @image_dir = $ftp->dir();
  } else {
    print "Error in images directory " . $emdbid . ". Message: " . $ftp->message;
  }
  my @map_dir = ();
  if ($ftp->cwd($map_cwd)) {
    @map_dir = $ftp->dir();
  } else {
    print "Error in map directory " . $emdbid . ". Message: " . $ftp->message;
  }
  if ($#xml_dir == 0 && $#image_dir >= 0 && $#map_dir >= 0) {
    my $image_file_index = -1;
    for my $current_image_index (0 .. $#image_dir) {
      my $base_image_filename =
          get_filename_from_dir_entry($image_dir[$current_image_index]);
      # There are a few names that people use more frequently.
      # If none of them are present then we'll just take the first image.
      if ($base_image_filename eq $emdbid ||
          $base_image_filename eq ("EMD-" . $emdbid) ||
          $base_image_filename eq ("emd_" . $emdbid) ||
          $base_image_filename eq ("emd-" . $emdbid)) { 
        $image_file_index = $current_image_index;
        last;
      }
    }
    if ($image_file_index == -1) {
      $image_file_index = 0;
    }
    print $output_file "$emdbid" .
          "\t" . get_date_from_dir_entry($map_dir[0]) .
          "\t" . get_date_from_dir_entry($xml_dir[0]) .
          "\t" . get_date_from_dir_entry($image_dir[$image_file_index]) .
          "\t" . $map_cwd . get_filename_from_dir_entry($map_dir[0]) .
          "\t" . $xml_cwd . get_filename_from_dir_entry($xml_dir[0]) .
          "\t" . $image_cwd .
                 get_filename_from_dir_entry($image_dir[$image_file_index]) .
          "\n";
  } else {
    print "Error reading files for " . $emdbid . ". XML found: " .
          ($#xml_dir + 1) . " Images found: " . ($#image_dir + 1) .
          " Maps found: " . ($#map_dir + 1) . "\n";
  }
  $ftp->quit();
}

print "Finished.\n";
print localtime() . "\n";


close($output_file);
