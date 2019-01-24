#include <ctime>
#include <sstream>
#include <iomanip>
#include <limits>
#include <fstream>
#include <ctime>

using namespace std;

#include <ANN/ANN.h>
#include "md.h"
#include "mdga_stats_options.h"
#include "mdgraph.h"
#include "pdb.h"
#include "transformations.h"
#include "statssubunits.h"

int main(int argc, char** argv)
{
	mdga_stats_options options(argc,argv);

	if(!options.parse_successful())
	{
		cerr << "Usage: multilzerd_subunit_stats --input <ga file> [-d <input_directory>] [--boundpath <directory>]  [-n number_of_predictions] [--timer]" << endl;
		exit(EXIT_FAILURE);
	}


	string ga_file = options.get_input_file();
	string directory = options.get_input_directory();
	string bound_directory = options.get_bound_directory();
	int results_number = options.get_results_number();
	decoy_program_t decoy_program_type;

	// Read the information on the file
	ifstream ga_file_stream(ga_file.c_str(), ifstream::in);
	string main_complex_name, chains_arg, discard;
	int generations, population_size;
	double clash_threshold;
	
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

	// we need to load both versions of the pdbs (with and without hydrogens)
	// because we only need the hydrogen-added version for scoring purposes
	vector<pdb> pdbs = load_pdbs(main_complex_name, chain_ids, directory, false);
	vector<pdb> bound_pdbs = load_pdbs(main_complex_name, chain_ids, bound_directory, false);
	vector<pdb> pdbs_with_hydrogen = load_pdbs(main_complex_name, chain_ids, directory, true);
	vector< vector<transformations*> > predictions = load_predictions(chain_ids, directory, decoy_program_type);

	// we provide the pdbs with hydrogens here because this method updates the internal score of each individual
	vector<md_edge_set> population = load_population(ga_file_stream);

	// close the ga file because we don't need it anymore
	ga_file_stream.close();

	// Used to measure the amount of time it takes to execute the program
	time_t start_time, end_time;
	start_time = time(NULL);

	///Get the contact information for the original complexes
	vector< vector<atom> > original_proteins;
	// Extract the atoms from the original pdb information
	for(size_t pdb_index = 0; pdb_index < bound_pdbs.size(); pdb_index++)
	{
		original_proteins.push_back(bound_pdbs[pdb_index].atoms);
	}

	// create the object that computes subunit statistics
	statssubunits stats_calculator(original_proteins, chain_ids);

	// we check of the number of results generated is limited
	size_t n = results_number == -1 ? population.size() : results_number;

	for(size_t population_index = 0; population_index < n; population_index++)
	{
		md_edge_set current_individual = population[population_index];

		// Get the transformed atoms according to the current individual's information
		vector< vector<atom> > decoy = current_individual.get_transformed_atoms(pdbs, predictions);
		
		vector<docking_stats> stats = stats_calculator.calculate_stats(decoy, population_index);
		// print stats for each subset
		for(size_t stats_index = 0; stats_index < stats.size(); stats_index++) {
			bool print_header = (!population_index) && (!stats_index); //i.e. if they are both zero
			cout << stats[stats_index].to_comma_string(print_header);
		}
	}

	end_time = time(NULL);

	if(options.get_use_timer()) {
		cerr << "Finished in " << (difftime(end_time, start_time) / 60.0) << " minutes" << endl;
	}

	
	//Explanation taken from the documentation
	//The ANN library allocates a small amount of storage, which is shared by all search structures
	//built during the program's lifetime. Because the data is shared, it is not deallocated,
	//even when the all the individual structures are deleted. To avoid the resulting (minor) memory
	//leak, the following function can be called after all search structures have been destroyed.
	annClose();
	exit(EXIT_SUCCESS);
}
