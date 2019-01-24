#!/usr/bin/perl

if ($#ARGV != 3 && $#ARGV != 4 ) {
  die "Usage: calculate_volume_and_3dzd.pl <contour level file> <std dev file>" .
      " <EM maps directory> <output volume file> [<output dir prefix>]\n\n" .
      "Example: calculate_volume_and_alternate_contours.pl contour_levels.txt std_dev.txt " .
      " map/ volume.txt [./] \n\nNote that 3DZD's are optionally generated if" .
      " the output directory is provided. Otherwise only volumes are calculated\n";
}

$contour_input_file = $ARGV[0];
$std_dev_file = $ARGV[1];
$maps_directory = $ARGV[2];
$volume_output_file = $ARGV[3];
if ($#ARGV == 4) {
  $compute_3dzd = 1;
  $output_dir_prefix = $ARGV[4];
} else {
  $compute_3dzd = 0;
}

# Clean up previously existing inv directories
if ($compute_3dzd) {
  `rm -Rf $output_dir_prefix/recommend_contour $output_dir_prefix/1_contour`;
  `rm -Rf $output_dir_prefix/2_contour $output_dir_prefix/1_2_contour`;
  `rm -Rf $output_dir_prefix/1_std_merge`;
  `mkdir $output_dir_prefix/recommend_contour; mkdir $output_dir_prefix/recommend_contour/inv`;
  `mkdir $output_dir_prefix/1_contour; mkdir $output_dir_prefix/1_contour/inv`;
  `mkdir $output_dir_prefix/2_contour; mkdir $output_dir_prefix/2_contour/inv`;
  `mkdir $output_dir_prefix/1_2_contour; mkdir $output_dir_prefix/1_2_contour/inv`;
  `mkdir $output_dir_prefix/1_std_merge; mkdir $output_dir_prefix/1_std_merge/inv`;
}

# Read all std dev values first and then retrieve them in the main loop that goes over all contour values/
%std_map = {};
open(STDDEVDATA, $std_dev_file);
foreach $stdline (<STDDEVDATA>) {
  @stdparts = split(/\s+/, $stdline);
  $std_map{$stdparts[0]} = $stdparts[1];
}
close(STDDEVDATA);

open(DATA, $contour_input_file);
open(VOLUME_FILE, ">$volume_output_file");

foreach $line (<DATA>) {
  @parts = split(/\s+/, $line);
  $emdbid = $parts[0];
  $contour = $parts[1];
  # default value that has the effect of generating the author recommended contour if it's missing
  # from the std dev input file.
  $std = 0;
  if (exists $std_map{$emdbid}) {
    $std = $std_map{$emdbid};
  }

  print "$emdbid\n";

  @em_command_results = `./em_volume $contour $maps_directory/$emdbid.map`;
  
  @volume_parts = split(/\s+/, $em_command_results[0]);
  $volume = $volume_parts[1];
  $max_density = substr($em_command_results[2],
                        rindex($em_command_results[2], " ") + 1);
  chomp($max_density);
  $one_third_contour = $contour + ($max_density - $contour) / 3;
  $two_thirds_contour = $contour + ($max_density - $contour) * (2 / 3);
  #$one_std_contour = $contour + $std;
  $one_std_contour = $std;
  # This generates the volume.txt file
  print VOLUME_FILE "$emdbid\t$volume\n";
  if ($compute_3dzd) {
    # Compute the three types of 3DZDs for different contour values
    `./map2zernike $maps_directory/$emdbid.map -n 20 -c $contour -p $emdbid-recommended`;
    `./map2zernike $maps_directory/$emdbid.map -n 20 -c $one_third_contour -p $emdbid-one-third`;
    `./map2zernike $maps_directory/$emdbid.map -n 20 -c $two_thirds_contour -p $emdbid-two-thirds`;
    `./map2zernike $maps_directory/$emdbid.map -n 20 -c $one_std_contour -p $emdbid-one-std`;
    # Make versions without the first line that contains the number of descriptors
    `tail -n +2 $emdbid-recommended.inv > $emdbid-recommended-raw.inv`;
    `tail -n +2 $emdbid-one-third.inv > $emdbid-one-third-raw.inv`;
    `tail -n +2 $emdbid-two-thirds.inv > $emdbid-two-thirds-raw.inv`;
    `tail -n +2 $emdbid-one-std.inv > $emdbid-one-std-raw.inv`;
    # Create the files for the 4 types of representations
    `cp $emdbid-recommended.inv $output_dir_prefix/recommend_contour/inv/$emdbid.inv`;
    `echo 242 > $output_dir_prefix/1_contour/inv/$emdbid.inv`;
    `cat $emdbid-recommended-raw.inv $emdbid-one-third-raw.inv >> $output_dir_prefix/1_contour/inv/$emdbid.inv`;
    `echo 242 > $output_dir_prefix/2_contour/inv/$emdbid.inv`;
    `cat $emdbid-recommended-raw.inv $emdbid-two-thirds-raw.inv >> $output_dir_prefix/2_contour/inv/$emdbid.inv`;
    `echo 363 > $output_dir_prefix/1_2_contour/inv/$emdbid.inv`;
    `cat $emdbid-recommended-raw.inv $emdbid-one-third-raw.inv $emdbid-two-thirds-raw.inv >> $output_dir_prefix/1_2_contour/inv/$emdbid.inv`;
    `echo 242 > $output_dir_prefix/1_std_merge/inv/$emdbid.inv`;
    `cat $emdbid-recommended-raw.inv $emdbid-one-std-raw.inv >> $output_dir_prefix/1_std_merge/inv/$emdbid.inv`;
    # Clean up all intermediate files for this EMDB ID
    `rm -f $emdbid-recommended* $emdbid-one-third* $emdbid-two-thirds* $emdbid-one-std*`;
  }
}
close(VOLUME_FILE);
close(DATA);
