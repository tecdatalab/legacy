#ifdef WITH_MPI

using namespace std;

#include "mpicoordinator.h"
#include "cluster_manager.h"

/*
 * Create a new mpicoordinator with the information provided
 */
mpicoordinator::mpicoordinator(vector<pdb>* proteins, vector< vector<transformations*> >* predictions,
						int clash_cutoff, double* weights)
{
	this->proteins = proteins;
	this->predictions = predictions;
	this->clash_cutoff = clash_cutoff;
	this->population_capacity = this->population_size = 0; //initially because we haven't received an update message yet
	this->weights = weights;
	// because the minimum spannning tree has vertices-1 edges
	this->num_edges = proteins->size() - 1;
	// allocate the memory used by the vector used for scoring purposes
	all_edges.reserve(num_edges);
	
	// memory used to receive messages
	edges_message = new mpi_message_update_score_send_t[num_edges];
	
	// setup the message types
	// create the mpi data type structures used for messaging
	this->main_to_aux_score = mpi_type_update_score_send();
	this->aux_to_main_score = mpi_type_update_score_receive();
	// a population update is basically a collection of edges
	this->main_to_aux_population_update = mpi_type_edge();
}

/*
 * Release the memory preallocated to receive messages
 */
mpicoordinator::~mpicoordinator()
{
	delete [] edges_message;
	
	if(this->population_capacity) { // then release the memory allocated
		delete [] population_message;
	}
	
	MPI_Type_free(&main_to_aux_score);
	MPI_Type_free(&aux_to_main_score);
	MPI_Type_free(&main_to_aux_population_update);
}

/*
 * When using MPI, this function will act as the main method for ranks
 * other than the main one (rank=0). This method will receive messages
 * from rank 0 containing the description of an edge set, it will calculate
 * clashes, physics and shape score, and it will send a message back to the main rank
 */
void mpicoordinator::distributed_process_main()
{
	MPI_Status status;

	do {
		// the message will always come from rank 0
		int error_code = MPI_Recv(edges_message, num_edges, main_to_aux_score, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if(error_code != MPI_SUCCESS) {
			cout << "Error detected while receiving general message on auxiliary process" << endl;
		}
		switch(status.MPI_TAG) {
			case MPI_MESSAGE_UPDATE_SCORE_SEND_TAG:
				process_update_score_message();
				break;
			case MPI_MESSAGE_UPDATE_PRE_POPULATION_TAG: // ignore message body and wait for a population update
				process_population_update_message();
				break;
			case MPI_MESSAGE_CLUSTER_SEND_TAG: //ignore the body too and proceed to cluster
				process_cluster_message();
				break;
		}
		
	} while(status.MPI_TAG != MPI_MESSAGE_FINISH_TAG);
}

/*
 * This method is called when the message received corresponds to
 * a score update request. It will do the appropriate calculations and send a message
 * back to rank 0 with the score values
 */
void mpicoordinator::process_update_score_message()
{
	// Will hold the information calculated for each of the decoys processed
	md_edge_set_score_t result;
	// This will have the same inherent values as result but it will have the format
	// specified by the MPI messaging header file
	mpi_message_update_score_receive_t result_msg;
	
	// reset the edges vector
	all_edges.clear();
	
	// transform the message into a collection of md_edge
	for(int msg_index = 0; msg_index < num_edges; msg_index++) {
		md_edge new_edge(edges_message[msg_index].receptor, edges_message[msg_index].ligand,
				edges_message[msg_index].prediction_number);
		all_edges.push_back(new_edge);
	}

	result = md_edge_set::update_score(&all_edges, *proteins, *predictions, clash_cutoff, weights);

	// update the values of the result message

	// extract the population index from the message
	result_msg.population_index = edges_message[0].population_index;

	result_msg.score = result.score;
	result_msg.pairwise_score = result.pairwise_score;
	result_msg.clashes = result.clashes;
	result_msg.soro_score = result.soro_score.get_base();

	// The message will be sent to rank 0
	MPI_Send(&result_msg, 1, aux_to_main_score, 0, MPI_MESSAGE_UPDATE_SCORE_RECEIVE_TAG, MPI_COMM_WORLD);
}

/*
 * After a pre-population update message was received, block
 * until the actual update comes and update this instance's population
 */
void mpicoordinator::process_population_update_message()
{
	MPI_Status status;
	int error_code, new_population_size;
	// first try to receive the size of the population
	error_code = MPI_Recv(&new_population_size, 1, MPI_INT, 0,
							MPI_MESSAGE_UPDATE_POPULATION_SIZE_TAG, MPI_COMM_WORLD, &status);
	if(error_code != MPI_SUCCESS) {
		cout << "Error detected while receiving population size message on auxiliary process" << endl;
	}
	
	// if the size increased the reallocate memory
	if(new_population_size > this->population_capacity) {
		if(this->population_capacity) {
			delete [] population_message;
		}
		this->population_message = new mpi_message_edge_t[this->num_edges * new_population_size];
		this->population_capacity = new_population_size;
	}
	// set this instances's population size
	this->population_size = new_population_size;

	// now that we know what size to expect, receive the actual info
	// each element has num_edges and we are receiving a total of population_size,
	// so the total received is num_edges*population_size
	error_code = MPI_Recv(this->population_message, this->num_edges * this->population_size,
								this->main_to_aux_population_update, 0,
								MPI_MESSAGE_UPDATE_POPULATION_TAG, MPI_COMM_WORLD, &status);
	if(error_code != MPI_SUCCESS) {
		cout << "Error detected while receiving population update message on auxiliary process" << endl;
	}
	
	// and proceed to update the instance's population
	// 1. for efficiency purposes check if the population vector has memory
	// assigned to it already. If it doesn't, reserve population_size 
	if(this->population.size() == 0) {
		this->population.reserve(this->population_size);
	}
	// 2. clear whatever we had before
	this->population.clear();
	// 3. create a local variable that will be used to create each individual before
	// adding it to the population
	md_edge_set current_individual;
	// 4. Each group of size num_edges in population_message will represent one individual
	// go through each of those groups and create an individual from them
	size_t message_index = 0;
	for(int population_index = 0; population_index < this->population_size; population_index++) {
		current_individual.clear();
		for(int edge_index = 0; edge_index < this->num_edges; edge_index++) {
			mpi_message_edge_t current_edge_info = this->population_message[message_index];
			md_edge new_edge(current_edge_info.receptor, current_edge_info.ligand,
							current_edge_info.prediction_number);
			current_individual.push_back(new_edge);
			message_index++;
		}
		this->population.push_back(current_individual);
	}
	/****** REMOVE this part 
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ofstream print_file;
	char numstring[10];
	sprintf(numstring, "%d", rank);
	string thename = string(numstring) + ".txt";
	print_file.open(thename.c_str());
	for(int population_index = 0; population_index < this->population_size; population_index++) {
		this->population[population_index].print(print_file);
	}
	print_file.close(); */
}

/*
 * Based on the rank of the auxiliary process, and the population
 * that should have been updated previously, calculated the RMSD's
 * corresponding to this process and send a message back to the main rank
 * with the results
 */
void mpicoordinator::process_cluster_message()
{
	int rank, numtasks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	
	// Calculate sizes for message passing purposes
	int auxranks = numtasks - 1;
	size_t n = this->population.size();
	size_t total_distances = (n * (n - 1)) / 2;
	size_t distances_per_rank = floor(total_distances / ((float)auxranks));
	
	//allocate a buffer that can be used to receive the distance updates from the auxiliary ranks
	double* distances = new double[distances_per_rank];
	
	size_t logical_start_index = distances_per_rank * (rank - 1);
	size_t logical_end_index = logical_start_index + (distances_per_rank - 1); //inclusive
	
	// create a cluster manager that uses a population instead of the explicit atomic structures
	cluster_manager cluster_handler(&(this->population));
	
	cluster_handler.calculate_distances(logical_start_index, logical_end_index, proteins, predictions, distances);
	
	// The message will be sent to rank 0
	MPI_Send(distances, distances_per_rank, MPI_DOUBLE, 0, MPI_MESSAGE_CLUSTER_RECEIVE_TAG, MPI_COMM_WORLD);
	
	delete [] distances;
}


/*
 *  This function will just send a termination message to all ranks different to 0
 *  The number of edges is used because, although we are not sending an edge set,
 *  the receiving end expects and edge set, so we will send an empty one
 */
void mpicoordinator::stop_distributed_processes(int num_edges)
{
	int numtasks;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	MPI_Datatype main_to_aux_type = mpi_type_update_score_send();
	mpi_message_update_score_send_t* edges_message = new mpi_message_update_score_send_t[num_edges];

	while(--numtasks > 0) {
		MPI_Send(edges_message, num_edges, main_to_aux_type, numtasks, MPI_MESSAGE_FINISH_TAG, MPI_COMM_WORLD);
	}

	MPI_Type_free(&main_to_aux_type);
	delete [] edges_message;
}

/*
 * This function sends a message to the rank identified by "process" in order to
 * calculate the score of the edge set provided
 *
 * The edges message pointer is provided to avoid allocating the buffer for the message
 * multiple times. edges_message points to an already allocated area that
 * can be used to copy the information in edges
 */
void mpicoordinator::send_distributed_score_request(md_edge_set edges, mpi_message_update_score_send_t* edges_message,
							MPI_Datatype main_to_aux_type, int population_index, int process)
{
	// copy the edges to the message structure used
	for(size_t edge_index = 0; edge_index < edges.size(); edge_index++) {
		md_edge current_edge = edges[edge_index];
		edges_message[edge_index].population_index = population_index;
		edges_message[edge_index].receptor = current_edge.receptor;
		edges_message[edge_index].ligand = current_edge.ligand;
		edges_message[edge_index].prediction_number = current_edge.prediction_number;
	}

	MPI_Send(edges_message, edges.size(), main_to_aux_type, process, MPI_MESSAGE_UPDATE_SCORE_SEND_TAG, MPI_COMM_WORLD);
}

/*
 * From the main rank, send messages to notify the auxiliary ranks that the population
 * has changed
 */
void mpicoordinator::send_distributed_population_update_request(vector<md_edge_set>* population)
{
	if(population->size()) {
		// Allocate and get the mpi types
		MPI_Datatype main_to_aux_type = mpi_type_update_score_send();
		MPI_Datatype population_update_type = mpi_type_edge();
	
		int population_size = population->size();
		size_t num_edges = population->at(0).size();
		mpi_message_edge_t* population_message = new mpi_message_edge_t[num_edges * population_size];
		// Each group of size num_edges in population_message will represent one individual
		// go through each of the individuals and create that consecutive-memory representation
		size_t message_index = 0;
		for(int population_index = 0; population_index < population_size; population_index++) {
			for(size_t edge_index = 0; edge_index < num_edges; edge_index++) {
				mpi_message_edge_t new_edge;
				md_edge original_edge = population->at(population_index)[edge_index];
				
				new_edge.receptor = original_edge.receptor;
				new_edge.ligand = original_edge.ligand;
				new_edge.prediction_number = original_edge.prediction_number;
				
				population_message[message_index++] = new_edge;
			}
		}
		
		// send a message with the size first and then a second message with the actual values
		// to each of the auxiliary ranks
		int numtasks;
		MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
		
		// used for the empty message next
		mpi_message_update_score_send_t* edges_message = new mpi_message_update_score_send_t[num_edges];
		
		for(int rank = 1; rank < numtasks; rank++) {
			// first, send an empty message indicating that the update is coming
			MPI_Send(edges_message, num_edges, main_to_aux_type, rank,
						MPI_MESSAGE_UPDATE_PRE_POPULATION_TAG, MPI_COMM_WORLD);
			// and then the actual info
			MPI_Send(&population_size, 1, MPI_INT, rank,
						MPI_MESSAGE_UPDATE_POPULATION_SIZE_TAG, MPI_COMM_WORLD);
			MPI_Send(population_message, num_edges * population_size, population_update_type, rank,
						MPI_MESSAGE_UPDATE_POPULATION_TAG, MPI_COMM_WORLD);
		}
	
		MPI_Type_free(&main_to_aux_type);
		MPI_Type_free(&population_update_type);
		delete [] edges_message;
		delete [] population_message;
	}
}

/*
 * From the main rank, send messages to notify the auxiliary ranks to start computing
 * clustering distances. The ranks know what to compute based on their rank id and
 * the population size
 * The number of edges is used because, although we are not sending an edge set,
 * the receiving end expects and edge set, so we will send an empty one
 */
void mpicoordinator::send_distributed_clustering_request(int num_edges)
{
	int numtasks;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	MPI_Datatype main_to_aux_type = mpi_type_update_score_send();
	mpi_message_update_score_send_t* edges_message = new mpi_message_update_score_send_t[num_edges];

	while(--numtasks > 0) {
		MPI_Send(edges_message, num_edges, main_to_aux_type, numtasks, MPI_MESSAGE_CLUSTER_SEND_TAG, MPI_COMM_WORLD);
	}

	MPI_Type_free(&main_to_aux_type);
	delete [] edges_message;
}

#endif
