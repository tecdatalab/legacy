#!/usr/bin/perl

# Assume that each of the parameters specified is a file and that
# we need to add an index to the beginning of each line, starting from zero

for $fileindex (0 .. $#ARGV)
{
	$currentfile = $ARGV[$fileindex];
	$rowindex = 0;
	
	# we'll keep the complete file in memory for I/O efficiency purposes
	# if this becomes a problem the program needs to be modified
	$newfiletext = "";

	
	for $currentline(`cat $currentfile`)
	{
		$newfiletext .= "$rowindex\t$currentline";
		$rowindex++;
	}

	#overwrite the file with the new text
	open(NEWFILE, ">$currentfile");

	print NEWFILE $newfiletext;

	close(NEWFILE);
}
