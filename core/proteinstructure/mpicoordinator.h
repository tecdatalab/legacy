#ifndef _MPICOORDINATOR_H_
#define _MPICOORDINATOR_H_

using namespace std;

#ifdef WITH_MPI

#include "mpimessaging.h"
#include "mdgraph.h"

/*
 * Instances of this class will be used to coordinate the
 * distributed process that perform computations requested by
 * the main MPI rank
 */
class mpicoordinator
{
	private:
		vector<pdb>* proteins;
		vector< vector<transformations*> >* predictions;
		int clash_cutoff;
		int population_size, population_capacity;
		double* weights;
		MPI_Datatype main_to_aux_score;
		MPI_Datatype aux_to_main_score;
		MPI_Datatype main_to_aux_population_update;
		int num_edges;
		// this vector will be the one passed to the scoring method. It will contain the edges
		// passed as the message, after transforming them to conform to the expected type
		vector<md_edge> all_edges;
		// stores the population information that may be updated using messages coming
		// from the main rank
		vector<md_edge_set> population;
		// This holds memory used to receive score update messages from the main rank
		mpi_message_update_score_send_t* edges_message;
		// Buffer that will store population update information
		mpi_message_edge_t* population_message;
	public:
		/*
		 * Create a new mpicoordinator with the information provided
		 */
		mpicoordinator(vector<pdb>* proteins, vector< vector<transformations*> >* predictions,
						int clash_cutoff, double* weights);
		
		/*
		 * Release the memory preallocated to receive messages
		 */
		 ~mpicoordinator();
		/*
		 * When using MPI, this function will act as the main method for ranks
		 * other than the main one (rank=0). This method will receive messages
		 * from rank 0 containing the description of an edge set, it will calculate
		 * clashes, physics and shape score, and it will send a message back to the main rank
		 */
		void distributed_process_main();
		
		/*
		 * This method is called when the message received corresponds to
		 * a score update request. It will do the appropriate calculations and send a message
		 * back to rank 0 with the score values
		 */
		void process_update_score_message();
		
		/*
		 * After a pre-population update message was received, block
		 * until the actual update comes and update this instance's population
		 */
		void process_population_update_message();
		
		/*
		 * Based on the rank of the auxiliary process, and the population
		 * that should have been updated previously, calculated the RMSD's
		 * corresponding to this process and send a message back to the main rank
		 * with the results
		 */
		void process_cluster_message();
									
		/*
		 * This function will just send a termination message to all ranks different to 0
		 * The number of edges is used because, although we are not sending an edge set,
		 * the receiving end expects and edge set, so we will send an empty one
		 */
		static void stop_distributed_processes(int num_edges);
		
		/*
		 * This function sends a message to the rank identified by "process" in order to
		 * calculate the score of the edge set provided
		 *
		 * The edges message pointer is provided to avoid allocating the buffer for the message
		 * multiple times. edges_message points to an already allocated area that
		 * can be used to copy the information in edges
		 */
		static void send_distributed_score_request(md_edge_set edges, mpi_message_update_score_send_t* edges_message,
									MPI_Datatype main_to_aux_type, int population_index, int process);
		/*
		 * From the main rank, send messages to notify the auxiliary ranks that the population
		 * has changed
		 */
		static void send_distributed_population_update_request(vector<md_edge_set>* population);
		/*
		 * From the main rank, send messages to notify the auxiliary ranks to start computing
		 * clustering distances. The ranks know what to compute based on their rank id and
		 * the population size
		 * The number of edges is used because, although we are not sending an edge set,
		 * the receiving end expects and edge set, so we will send an empty one
		 */
		static void send_distributed_clustering_request(int num_edges);

};
#endif

#endif
