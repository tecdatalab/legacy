#!/usr/bin/perl

# Receives a file that contains correlations between patches
# and adds the fractions of those patches that correspond to the interface.
# Since the correlations files contents are of the form: 
# 
# patchfile1.pdb.ms.inv patchfile2.pdb.ms.inv <correlation value>
#
# it is assumed that the current directory contains files
# patchfile1.pdb.ms.int patchfile2.pdb.ms.int that have all the atoms in the
# real interface.
#
# The original contents of the pdb.ms patches are output with 4
# additional columns:
#	1. number_of_patch_atoms_on_interface / total_patch_atoms (for patch1)
#	2. number_of_patch_atoms_on_interface / total_interface_atoms (for patch1)
#	3. Same as #1 for patch2
#	4. Same as #2 for patch2

if ($#ARGV + 1 != 1) {
	print "Usage: join_correlations_interface.pl <correlations file>\n";
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
		$patch1_int_file = "$patch1_file.int";
		$patch1_int_file =~ s/-atom.....residue....//;
		$patch2_int_file = "$patch2_file.int";
		$patch2_int_file =~ s/-atom.....residue....//;

		$patch1_interface_data =
				`../interface_in_patch.py $patch1_file $patch1_int_file`;
		$patch2_interface_data =
				`../interface_in_patch.py $patch2_file $patch2_int_file`;
		chomp($patch1_interface_data);
		chomp($patch2_interface_data);

		print "$corr_line $patch1_interface_data $patch2_interface_data\n";
	}
}
