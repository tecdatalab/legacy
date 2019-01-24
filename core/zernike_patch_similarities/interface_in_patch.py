#!/usr/bin/python

import sys

import atom

def main():
	"""Receives a PDB with a patch and a PDB that contains only the interface
	and computes two numbers:
	1. number_of_patch_atoms_on_interface / total_patch_atoms
	2. number_of_patch_atoms_on_interface / total_interface_atoms
	Given a mark_sur-modified PDB file, generate one patch per surface
	residue that includes all atoms closer than the threshold given.
	A new PDB file is written for each patch.
	"""
	if len(sys.argv) == 3:

		patch_pdb = atom.AtomCollection(sys.argv[1], atom.Atom.INPUT_PDB_BASIC)
		interface_pdb = atom.AtomCollection(sys.argv[2], atom.Atom.INPUT_PDB_BASIC)

		overlap_count = float(interface_pdb.count_overlap(patch_pdb))

		fraction_wrt_patch_total = overlap_count / patch_pdb.length()
		fraction_wrt_interface_total = overlap_count / interface_pdb.length()

		print str(fraction_wrt_patch_total) + " " + \
					str(fraction_wrt_interface_total)

	else:
		print "Usage: " + sys.argv[0] + " <patch pdb> <interface pdb>"

# Just call the main function.
main()

