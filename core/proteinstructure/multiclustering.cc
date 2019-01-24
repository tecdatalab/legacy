#include <iostream>
#include <fstream>
#include <vector>
#include "multiclustering_options.h"
#include "mdstate.h"

#ifdef WITH_MPI

// only include if the mpi version is being compiled
#include "mpi.h"
#include "mpicoordinator.h"

#endif

using namespace std;

int main(int argc, char** argv)
{
#ifdef WITH_MPI
	
	/*** Perform the MPI Initialization that applies to all ranks ***/
	int rank, numtasks;
	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

#endif

	multiclustering_options options(argc,argv);
	if(!options.parse_successful())
	{
		cerr << "Usage: multilzerd_cluster --input <MultiLZerD output file> -o <output filename> --cutoff <cluster distance>"
			<< " [--printsizes <cluster sizes output file>] [-d <input_directory>]" << endl;
		exit(EXIT_FAILURE);
	}

	time_t start_time, end_time;
	start_time = time(NULL);

	ofstream output_file;
	ofstream sizes_file;

	string input_filename = options.get_input_file();
	string output_filename = options.get_output_file();
	string directory = options.get_input_directory();
	string sizes_filename = options.get_sizes_file();
	bool print_sizes = !sizes_filename.empty();
	double cutoff = options.get_cutoff();

	mdstate population = mdstate(input_filename, directory, cutoff);

#ifdef WITH_MPI
	/*** This is the MPI code that will be executed for ranks other than 0 ***/
	if(rank != 0) {
		// This will just call the auxiliary procedures and after they finish they will exit
		// Note: zero means there is no clash cutoff that aborts the counting process
		mpicoordinator coordinator(population.get_pdbs(),
									population.get_predictions(),
									population.get_clash_cutoff(),
									population.get_weights());
		coordinator.distributed_process_main();

		// Wrap up MPI things
		MPI_Finalize();

		return 0;
	}

#endif

	vector<cluster_neighborhood> cluster_centers = population.cluster();

	/* write the output after clustering */

	output_file.open(output_filename.c_str());
	population.write(output_file);
	output_file.close();

	if(print_sizes)
	{
		// delete it first, then open it
		remove(sizes_filename.c_str());
		sizes_file.open(sizes_filename.c_str());
	}

	for(size_t center_index = 0; center_index < cluster_centers.size(); center_index++)
	{
		sizes_file << cluster_centers[center_index].get_neighbor_count() << endl;
	}

	// close the sizes file in case we used it
	if(print_sizes)
	{
		sizes_file.close();
	}

	end_time = time(NULL);

	cout << "Finished in " << (difftime(end_time, start_time) / 60.0) << " minutes" << endl;

#ifdef WITH_MPI

	// Wrap up MPI things (only rank == 0 should get here)
	
	mpicoordinator::stop_distributed_processes(population.get_pdbs()->size() - 1);

	MPI_Finalize();

#endif

	exit(EXIT_SUCCESS);
}
