#!/usr/bin/python

import math
import string
import sys

import atom

def main():
	"""Given two PDB files, return the distance between the centroids
	of each PDB file. This takes into account ALL ATOM elements in the
	PDB.
	"""
	if len(sys.argv) == 3:
		filename1 = sys.argv[1]
		filename2 = sys.argv[2]

		collection1 = atom.AtomCollection(filename1, atom.Atom.INPUT_PDB_BASIC)
		collection2 = atom.AtomCollection(filename2, atom.Atom.INPUT_PDB_BASIC)
		result = str(collection1.calculate_centroid_distance(collection2))
		
		# This section is targeted towards the patch experiment
		# done for a grant in July 2012 where the filenames contained
		# the number of the atom that was used as patch center
		atom1_index = string.find(filename1, "atom") + 4
		atom2_index = string.find(filename1, "atom") + 4

		if atom1_index != -1 and atom2_index != -1:
			serial1 = int(filename1[atom1_index:atom1_index + 5])
			serial2 = int(filename2[atom2_index:atom1_index + 5])

			atom1 = collection1.get_atom_by_serial_number(int(serial1))
			atom2 = collection2.get_atom_by_serial_number(int(serial2))

			deltax = atom1.x - atom2.x
			deltay = atom1.y - atom2.y
			deltaz = atom1.z - atom2.z
			
			result += " " + str(math.sqrt(deltax * deltax +
																		deltay * deltay +
																		deltaz * deltaz))

		print result

	else:
		print "Usage: " + sys.argv[0] + " <pdb1> <pdb2>"

# Just call the main function.
main()

