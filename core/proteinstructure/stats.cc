#include <ctime>
#include <sstream>
#include <iomanip>
#include <limits>
#include <fstream>
#include <ctime>

using namespace std;

#include "align.h"
#include "contact.h"
#include "md.h"
#include "mdga_stats_options.h"
#include "mdgraph.h"
#include "rmsd.h"
#include "pdb.h"
#include "transformations.h"
#include "utils.h"

#ifdef WITH_MPI

// only include if the mpi version is being compiled
#include "mpi.h"
#include "mpicoordinator.h"

#endif

int main(int argc, char** argv)
{
#ifdef WITH_MPI
	
	/*** Perform the MPI Initialization that applies to all ranks ***/
	int rank, numtasks;
	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

#endif
	mdga_stats_options options(argc,argv);

	if(!options.parse_successful())
	{
		cerr << "Usage: multilzerd_stats --input <ga file> [-d <input_directory>] [--boundpath <directory>] "
			<< "[-n number_of_predictions] [--timer] [--noscore] [--calphaonly]"
			<< WEIGHT_PARAM_DESCRIPTION << endl;
		exit(EXIT_FAILURE);
	}


	string ga_file = options.get_input_file();
	string directory = options.get_input_directory();
	string bound_directory = options.get_bound_directory();
	int results_number = options.get_results_number();
	decoy_program_t decoy_program_type;
	double* weights = options.get_weights();

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
	
	vector<pdb> pdbs = load_pdbs(main_complex_name, chain_ids, directory, false, options.get_calpha_only());
	vector<pdb> bound_pdbs = load_pdbs(main_complex_name, chain_ids, bound_directory, false, options.get_calpha_only());
	vector<pdb> pdbs_with_hydrogen;
	// initialize only if scores are generated
	if(options.get_generate_scores()) {
		pdbs_with_hydrogen = load_pdbs(main_complex_name, chain_ids, directory, true);
	}
	vector< vector<transformations*> > predictions = load_predictions(chain_ids, directory, decoy_program_type);

	// we provide the pdbs with hydrogens here because this method updates the internal score of each individual
	vector<md_edge_set> population = load_population(ga_file_stream);

	// close the ga file because we don't need it anymore
	ga_file_stream.close();

	// Used to measure the amount of time it takes to execute the program
	time_t start_time, end_time;
	start_time = time(NULL);

#ifdef WITH_MPI
	
	/*** This is the MPI code that will be executed for ranks other than 0 ***/
	if(rank != 0) {
		// This will just call the auxiliary procedures and after they finish they will exit
		// Note: zero means there is no clash cutoff that aborts the counting process

		// Auxiliary processes only score, thus if no scoring is performed just call MPI_Finalize

		if(options.get_generate_scores()) {
			mpicoordinator coordinator(&pdbs_with_hydrogen, &predictions, 0 , weights);
			coordinator.distributed_process_main();
		}

		// Wrap up MPI things
		MPI_Finalize();

		annClose();
		exit(EXIT_SUCCESS);
	}

#endif
	///Get the contact information for the original complexes
	vector< vector<atom> > original_proteins;
	// Extract the atoms from the original pdb information
	for(size_t pdb_index = 0; pdb_index < bound_pdbs.size(); pdb_index++)
	{
		original_proteins.push_back(bound_pdbs[pdb_index].atoms);
	}
	// and now get the contact information
	vector<string> original_contacts;
	vector< vector<string> > original_contacts_per_chain;
	identify_contact_residues(original_proteins, original_contacts, original_contacts_per_chain);
	//...as well as marking the interface residues
	get_interface_residues(original_proteins);

	//This array will hold the rmsd values for each chain
	double* rmsd_by_chain = new double[chain_ids.size()];

	// Write the column headers
	cout << "rank,fnat,fnon-nat";
	for(size_t chain_index = 0; chain_index < chain_ids.size(); chain_index++)
	{
		cout << ",fnat " << chain_ids[chain_index] << ",fnon-nat " << chain_ids[chain_index];
	}
	cout << ",Avg fnat,Avg fnon-nat,RMSD";
	for(size_t chain_index = 0; chain_index < chain_ids.size(); chain_index++)
	{
		cout << ",RMSD " << chain_ids[chain_index];
	}
	cout << ",Avg RMSD,iRMSD";
	for(size_t chain_index = 0; chain_index < chain_ids.size(); chain_index++)
	{
		cout << ",iRMSD " << chain_ids[chain_index];
	}
	cout << ",Avg iRMSD,PairwiseScore,Clashes";
	if(options.get_generate_scores()) {
		cout << ",Score,vdw,vdw_attr,vdw_rep,elec,elec_sr_attr,elec_lr_attr,elec_sr_rep,elec_lr_rep,hbp_ss,solv,sasa,acp";
	}
	cout << endl;
	// average variables
	double avg_fnat, avg_fnon_nat, avg_rmsd, avg_irmsd;
	// variable used for fnat storage, they are reset on each iteration
	double fnat, fnon_nat;
	double* fnat_per_chain = new double[chain_ids.size()];
	double* fnon_nat_per_chain = new double[chain_ids.size()];

	// we check of the number of results generated is limited
	size_t n = results_number == -1 ? population.size() : results_number;


#ifdef WITH_MPI
	/* Update the scores in the population first. This method will take advantage of
	 * mpi parallelization if available
	 * the zero represents the clash cutoff; when zero is provided as cutoff
	 * it is actually ignored, which is correct in this case because we only want to calculate
	 * the number of clashes without trying to discriminate between different conformations
	 */
	if(options.get_generate_scores()) {
		md_edge_set::update_score(&population, pdbs_with_hydrogen, predictions, 0, weights, n);
	}
#endif
	for(size_t population_index = 0; population_index < n; population_index++)
	{
		// reset the averages to zero
		avg_fnat = avg_fnon_nat = avg_rmsd = avg_irmsd = 0;
		// This map is used to calculate the RMSD
//		map<pair<int,int>, vector<int> > M;
		vector<string> individual_contacts;
		vector< vector<string> > individual_contacts_per_chain;
		
		md_edge_set current_individual = population[population_index];

#ifndef WITH_MPI
		// if we are not using mpi calculate it directly
		if(options.get_generate_scores()) {
			current_individual.update_score(pdbs_with_hydrogen, predictions, 0, weights);
		}
#endif
		
		// Get the transformed atoms according to the current individual's information
		vector< vector<atom> > transformed_proteins = current_individual.get_transformed_atoms(pdbs, predictions);
		
		// calculate fnat
		identify_contact_residues(transformed_proteins, individual_contacts, individual_contacts_per_chain);
		get_fnat(original_contacts, individual_contacts, original_contacts_per_chain, individual_contacts_per_chain,
				fnat, fnon_nat, fnat_per_chain, fnon_nat_per_chain);
		
		//...and then the RMSD
		double rmsd = calculate_multiple_rmsd(original_proteins, transformed_proteins, rmsd_by_chain, false);
//		int n = get_matching_atoms(transformed_proteins, original_proteins, M);
//		double rmsd = align(transformed_proteins, original_proteins, M, n);

		// first, print the rank of the individual in the population
		cout << (population_index + 1);

		// print the metrics for each individual
		cout << "," << fnat << "," << fnon_nat;
		// print the fnat/fnon-nat for each chain
		for(size_t fnat_index = 0; fnat_index < chain_ids.size(); fnat_index++)
		{
			// update averages
			avg_fnat += fnat_per_chain[fnat_index];
			avg_fnon_nat += fnon_nat_per_chain[fnat_index];

			cout << "," << fnat_per_chain[fnat_index] << "," << fnon_nat_per_chain[fnat_index];
		}

		// divide by the number of chains
		avg_fnat /= chain_ids.size();
		avg_fnon_nat /= chain_ids.size();
		cout << "," << avg_fnat << "," << avg_fnon_nat;
		
		cout << "," << rmsd;
		// print the rmsd for each chain
		for(size_t rmsd_index = 0; rmsd_index < chain_ids.size(); rmsd_index++)
		{
			// update average
			avg_rmsd += rmsd_by_chain[rmsd_index];
			cout << "," << rmsd_by_chain[rmsd_index];
		}

		//divide by the number of chains
		avg_rmsd /= chain_ids.size();
		cout << "," << avg_rmsd;

		// as well as the irmsd
		double irmsd = calculate_multiple_rmsd(original_proteins, transformed_proteins, rmsd_by_chain, true);
		cout << "," << irmsd;
		for(size_t rmsd_index = 0; rmsd_index < chain_ids.size(); rmsd_index++)
		{
			// update average
			avg_irmsd += rmsd_by_chain[rmsd_index];
			cout << "," << rmsd_by_chain[rmsd_index];
		}
		//divide by the number of chains
		avg_irmsd /= chain_ids.size();

		cout << "," << avg_irmsd << "," << current_individual.get_pairwise_score()
			<< "," << current_individual.get_clashes();
		if(options.get_generate_scores()) {
			cout << "," << current_individual.get_score() << ",";
			current_individual.print_detailed_score(cout);
		}
		cout << endl;
	}

	end_time = time(NULL);

	if(options.get_use_timer()) {
		cerr << "Finished in " << (difftime(end_time, start_time) / 60.0) << " minutes" << endl;
	}

	delete[] rmsd_by_chain;
	delete[] fnat_per_chain;
	delete[] fnon_nat_per_chain;
	
	

#ifdef WITH_MPI

	// Wrap up MPI things (only rank == 0 should get here)
	// Processes were triggered only if scoring was enabled
	if(options.get_generate_scores()) {
		mpicoordinator::stop_distributed_processes(pdbs.size() - 1);
	}

	MPI_Finalize();

#endif
	//Explanation taken from the documentation
	//The ANN library allocates a small amount of storage, which is shared by all search structures
	//built during the program's lifetime. Because the data is shared, it is not deallocated,
	//even when the all the individual structures are deleted. To avoid the resulting (minor) memory
	//leak, the following function can be called after all search structures have been destroyed.
	annClose();
	exit(EXIT_SUCCESS);
}
