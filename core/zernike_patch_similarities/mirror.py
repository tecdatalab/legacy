#!/usr/bin/python

import sys

import atom

def main():
  """Given a PDB file, generate two new files, one shifted to the origin
  while the second one will be the mirror image flipping x-coordinates.
  """
  if len(sys.argv) == 4:
    input_file = sys.argv[1]
    origin_file = sys.argv[2]
    mirror_file = sys.argv[3]

    input_pdb = atom.AtomCollection(input_file, atom.Atom.INPUT_PDB_BASIC)

    # Returns a copy of the atoms after shifting them to the origin
    origin_pdb = input_pdb.get_atoms_shifted_to_origin()
    # Same here, a copy after flipping x-coordinates
    mirror_pdb = origin_pdb.get_x_mirrored_atoms()

    origin_pdb.write_to_file(origin_file)
    mirror_pdb.write_to_file(mirror_file)
    
  else:
    print "Usage: " + sys.argv[0] + " <input pdb> <origin-shifted pdb> " +\
          "<mirror pdb>"

# Just call the main function.
main()

