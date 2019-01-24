#!/usr/bin/python

import sys

import atom

def main():
	"""Given a series of PDB files, it will create a new file that contains
	the atoms on the interface between every pairwise combination of them.
	Each file created will have the original filename with a ".int" suffix.
	"""
	if len(sys.argv) >= 3:
		filenames = sys.argv[1:]
		interfaces = \
				atom.AtomCollection.calculate_interfaces_from_files(\
						filenames, atom.Atom.INPUT_PDB_BASIC)
		for file_index in range(len(filenames)):
			interfaces[file_index].write_to_file(filenames[file_index] + ".int")
	else:
		print "Usage: " + sys.argv[0] + " <pdb1> <pdb2> ... <pdbn>"

# Just call the main function.
main()

