#!/usr/bin/python

import sys

import atom

def main():
	"""Given a mark_sur-modified PDB file, generate one patch per surface
	residue that includes all atoms closer than the threshold given.
	A new PDB file is written for each patch.
	"""
	if len(sys.argv) == 3:
		filename = sys.argv[1]
		threshold = int(sys.argv[2])
		collection = atom.AtomCollection(filename, atom.Atom.INPUT_MARK_SUR)
		patches = collection.get_surface_patches(threshold)
		for atom_collection in patches:
			atom_collection.write_to_file()
	else:
		print "Usage: " + sys.argv[0] + " <mark_sur pdb> <distance threshold>"

# Just call the main function.
main()

