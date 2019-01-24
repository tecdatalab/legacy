#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include "pairwisestats_options.h"
#include "contact.h"
#include "lzerd_transformations.h"
#include "pdb.h"
#include "rmsd.h"

using std::cerr;
using std::endl;
using std::string;

int main(int argc, char** argv)
{
	pairwisestats_options options(argc,argv);
	if(!options.parse_successful())
	{
		cerr << "Usage: lzerd_stats -R <receptor pdb file> -L <ligand pdb file> --pdb <predicted ligand in PDB format>" << endl;
		exit(EXIT_FAILURE);
	}

	string pdb_file = options.get_pdb_file();
	string ligand_file = options.get_ligand_file();
	string receptor_file = options.get_receptor_file();

	// Load the information
	pdb receptor, ligand, prediction;

	read_protein(ligand_file, ligand);
	read_protein(receptor_file, receptor);
	read_protein(pdb_file, prediction);

	// Data structures needed to calculate LRMSD
	vector<string> LRES;
	get_residue_list(ligand.atoms, LRES);

	int k3 = 0;
	map<int, int> MRES_UL; //matching residues
	set_ligand_CA(MRES_UL, LRES, ligand.atoms, prediction.atoms, k3);
	check_ligand_CA(MRES_UL);

	double lrmsd = calculate_lrms(ligand.atoms, prediction.atoms, MRES_UL);

	cout << setprecision(3) << left << setw(10) << lrmsd << endl;

	exit(EXIT_SUCCESS);
}
