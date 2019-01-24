#!/usr/bin/perl
if($#ARGV + 1 != 4 && $#ARGV +1 != 5)
{
	print "This program takes a situs volume file, expands the data into a cube (by adding zeroes) and replaces every value by a 0 or a 1\n";
	print "If a voxel value is within the range provided, it's replaced by a 1 or a 0 otherwise\n";
	print "The main purpose of this program is to generate a repreentation of the file that can be used by the zernike program (in \"grid\" format)\n";
	print "An optional 4th parameter can be supplied in order to generate a situs file instead of a grid (for visualization purposes)\n\n";
	print "Usage: situs2grid.pl <input.situs file> <output.grid file> <lower bound> <upper bound> [--situs]\n";
}
else
{
	$input_path = $ARGV[0];
	$output_path = $ARGV[1];
	$lower_bound = $ARGV[2];
	$upper_bound = $ARGV[3];
	$use_situs_format = $ARGV[4] eq "--situs";
	$floating_format = "%11.6f";
	open(INPUT, $input_path);
	open(OUTPUT, ">$output_path");

	# read the header of the situs file
	$header_line = <INPUT>;
	chomp($header_line);
	@split_header = split(/ /, $header_line);
	# this is the resolution of the map in angstroms
	$voxel_width = $split_header[0];
	# x,y,z coordinates of the first voxel, it also represents the origin of the map
	$originx = $split_header[1];
	$originy = $split_header[2];
	$originz = $split_header[3];
	# number of voxels on each axis, the total of voxels is given by nx * ny * nz
	$nx = $split_header[4];
	$ny = $split_header[5];
	$nz = $split_header[6];
	# we need to convert the situs file to a cubic grid (all dimensions are the same size)
	# we get the max of the 3 dimensions to fill the rest of the cube
	$max = $nx < $ny ? $ny : $nx;
	$max = $nz < $max ? $max : $nz;

	#ignore the first line since its blank
	$next_line = <INPUT>;

	if($use_situs_format) # print the header, otherwise don't because the grid format has only numbers
				# and it assumes that the box is a cube
	{
		print(OUTPUT "$voxel_width $originx $originy $originz $max $max $max\n\n");
	}

	$voxel_number = 1;
	# we want to have a cube with the original values so we will just fill up a
	# cube by adding zeroes to the two "shorter" dimensions
	for($z = 0; $z < $max; $z++)
	{
		for($y = 0; $y < $max; $y++)
		{
			for($x = 0; $x < $max; $x++)
			{
				# check if the indices are outside the original range, in which case we need to print a 0
				if($z >= $nz || $y >= $ny || $x >= $nx)
				{
					$voxel_value = 0;
				}
				else # it's within the range so take the next value from the file
				{
					if($#elements < 0) # meaning that the length is zero, we need to load more from the file
					{
						$next_line = <INPUT>;
						chomp($next_line);
						@elements = split(/ +/, $next_line);
						shift(@elements); # the first element after the split is a blank so we remove it
					}
					# remove the next element from the input stream (stored in @elements array)
					# and use that as the next voxel value
					$original_voxel_value = shift(@elements);
					# now determine if it's 0 or 1 based on the threshold
					$voxel_value = ($lower_bound < $original_voxel_value && $original_voxel_value < $upper_bound) ? 1 : 0;
				}
				if($use_situs_format) # then use the floating point formatting
				{
					printf(OUTPUT "$floating_format ", $voxel_value);
				}
				else # just print 0 or 1
				{
					print(OUTPUT "$voxel_value ");
				}

				if($voxel_number % 10 == 0 && $use_situs_format) #print a newline every 10 voxels (the same as the situs format)
								# only print newlines if using situs format
				{
					print(OUTPUT "\n");
				}
				$voxel_number++;
			}
		}
	}

	close(INPUT);
	close(OUTPUT);
}
