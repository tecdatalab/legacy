#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "lzerd_transformations.h"
#include "pdb.h"

using namespace std;

int main(int argc, char** argv)
{
  if(argc < 5) {
    cout << "Usage: PDBGEN receptor.pdb ligand.pdb prediction_file "
         << "number_of_predictions_to_output" << endl;
    exit(EXIT_FAILURE);
  }

  string receptor_filename = argv[1];
  string ligand_filename = argv[2];
  string lzerd_filename = argv[3];
  size_t predictions_output = static_cast<size_t>(atoi(argv[4]));

  pdb receptor, ligand;
  read_protein(receptor_filename, receptor);
  read_protein(ligand_filename, ligand);

  lzerd_transformations lzerd_predictions(lzerd_filename);

  vector<atom> transformed;
  char current_filename[20];

  for(size_t ligand_index = 0;
      ligand_index < lzerd_predictions.get_size()
          && ligand_index < predictions_output; ligand_index++) {
    transformed.clear();
    lzerd_predictions.transform_atoms(ligand.atoms, transformed, ligand_index);
    sprintf(current_filename, "ligand%d.pdb",
            static_cast<int>(ligand_index + 1));
    write_complex(current_filename, transformed);
  }

  exit(EXIT_SUCCESS);
}
