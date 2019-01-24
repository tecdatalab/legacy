#include <ctime>
#include <sstream>
#include <iomanip>
#include <limits>
#include <fstream>

using namespace std;

#include "md.h"
#include "mdgraph.h"
#include "pdb.h"
#include "transformations.h"
#include "utils.h"

int main(int argc, char** argv)
{
	if(argv[1] == NULL || argv[2] == NULL || argv[3] == NULL || argv[4] == NULL || argv[5] == NULL)
	{
		cerr << "Usage: multilzerd_create_pdb ga_file input_directory start end pdb_prefix" << endl;
		exit(EXIT_FAILURE);
	}

	string ga_file(argv[1]), directory(argv[2]), start_str(argv[3]), end_str(argv[4]), pdb_prefix(argv[5]);
	size_t start_index = atoi(start_str.c_str()) - 1; /* -1 because the input is not zero-indexed */
	size_t end_index = atoi(end_str.c_str()) - 1;

	// Read the information on the file
	ifstream ga_file_stream(ga_file.c_str(), ifstream::in);
	string main_complex_name, chains_arg, discard;
	int generations, population_size;
	double clash_threshold;
	decoy_program_t decoy_program_type;
	string decoy_program;
	string score_type;
	// do it twice because we have "PDBID 1VCB", for example. The first one read PDBID to discard it
	ga_file_stream >> discard;
	ga_file_stream >> main_complex_name;
	// twice for the same reason
	ga_file_stream >> discard;
	ga_file_stream >> chains_arg;
	vector<string> chain_ids = split(chains_arg, ',');
	// Number of generations
	ga_file_stream >> discard;
	ga_file_stream >> generations;
	// Population size
	ga_file_stream >> discard;
	ga_file_stream >> population_size;
	// Clash threshold
	ga_file_stream >> discard;
	ga_file_stream >> clash_threshold;
	// Score type
	ga_file_stream >> discard;
	ga_file_stream >> score_type;
	// Decoy program
	ga_file_stream >> discard;
	ga_file_stream >> decoy_program;

	decoy_program_type = decoy_program.compare(DECOYPROGRAM_LZERD) == 0 ?
	                                        decoy_program_lzerd : decoy_program_zdock;
	
	// The last parameter is set to false because when we create the PDB files we don't want
	// hydrogens added to them
	vector<pdb> pdbs = load_pdbs(main_complex_name, chain_ids, directory, false);
	vector< vector<transformations*> > predictions = load_predictions(chain_ids, directory, decoy_program_type);
	vector<md_edge_set> population = load_population(ga_file_stream);


	// close the ga file because we don't need it anymore
	ga_file_stream.close();

	// Go from start to end indices and write the pdb complex to the current directory
	for(;start_index <= end_index && start_index < population.size(); start_index++)
	{
		population[start_index].write(pdbs, predictions, pdb_prefix, start_index);
	}

	exit(EXIT_SUCCESS);
}
