#!/usr/bin/python

import sys

import atom

def main():
  """Analyze if the residues given are on the interface.
  It's assumed that a single PDB file contains two chains
  Params:
    PDB file
    Residue list: Follows a comma-separated CHAINLETTER-### format
      (e.g. A-55,A-58,B-66,B-101)
  """
  if len(sys.argv) >= 3:
    atom_collections = AtomCollection.new_collections_by_chain(sys.argv[1])
    if len(atom_collections) != 2:
      print "Two chains were expected in " + sys.argv[0]
      return None
    else:
      receptor = atom_collections[0]
      ligand = atom_collections[1]
      receptor.set_atoms_in_interface(ligand)
      ligand.set_atoms_in_interface(receptor)

      residue_list = sys.argv[2].split(",")

      print residue_list
      

  else:
    print "Usage: " + sys.argv[0] + \
          " <pdb> <residue list> (e.g. A-55,A-58,B-66,B-101)"

# Just call the main function.
main()

