#include "mdstate.h"
#include "mdgraph.h"
#include "mdga_options.h"


#include <fstream>
#include <ctime>

using std::cout;
using std::endl;

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

	ofstream log_file;
	ofstream output_file;

	mdga_options options(argc,argv);
	if(!options.parse_successful())
	{
		cout << "Usage: multilzerd --pdbid <PDB ID> --chains <A,B,C,etc> -o <output_prefix>" << endl
         << "\t[--interactions <chainX-chainY,chainH-chainM,...>] [--generations <repetitions>]" << endl
         << "\t[--population <size>] [--clashes <num>] [-d <input_directory>]" << endl
         << "\t[--cluster <rmsd threshold>] [--renew <new random elements per generation>] " << endl
         << "\t[--stagnation <num>] [--crossover] [--detailed] [--resume output.ga.out] " << endl
         << "\t" << WEIGHT_PARAM_DESCRIPTION << endl;
		return 0;
	}
	
	mdstate current_generation = mdstate(&options);
	
#ifdef WITH_MPI
	/*** This is the MPI code that will be executed for ranks other than 0 ***/
	if(rank != 0) {
		// This will just call the auxiliary procedures and after they finish they will exit
		// Note: zero means there is no clash cutoff that aborts the counting process
		mpicoordinator coordinator(current_generation.get_pdbs(),
									current_generation.get_predictions(),
									current_generation.get_clash_cutoff(),
									current_generation.get_weights());
		coordinator.distributed_process_main();

		// Wrap up MPI things
		MPI_Finalize();

		return 0;
	}

#endif

	string log_filename = options.get_log_filename();
	string output_filename = options.get_output_filename();

	// Used to measure the amount of time it takes to execute the program
	time_t start_time, end_time;
	start_time = time(NULL);

	log_file.open(log_filename.c_str());

	/* update the current generation to the next until we finish (i.e. advance_generation returns false) */
	while(current_generation.advance_generation(log_file))
	{
		// If we need to write the detailed representation we invoke the method in current_generation
		if(options.is_detailed())
		{
			char numstring[10];
			sprintf(numstring, "%05d-", current_generation.get_current_generation());
			string detailed_filename = string(numstring) + output_filename;
			
			output_file.open(detailed_filename.c_str());
			current_generation.write(output_file);
			output_file.close();
		}

		// Since the detailed output is optional we still update the main output file after each generation
		log_file << "Writing generation's best candidates..." << endl;

		output_file.open(output_filename.c_str());
		current_generation.write(output_file);
		output_file.close();

		// close and open the log file again to delete its contents
		log_file.close();
		log_file.open(log_filename.c_str());
	}


	end_time = time(NULL);

	// the log file should already be open
	log_file << "Finished in " << (difftime(end_time, start_time) / 60.0) << " minutes" << endl;
	log_file.close();
	
#ifdef WITH_MPI

	// Wrap up MPI things (only rank == 0 should get here)
	
	mpicoordinator::stop_distributed_processes(current_generation.get_pdbs()->size() - 1);

	MPI_Finalize();

#endif
	
	return 0;
}
