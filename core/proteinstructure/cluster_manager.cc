#include "cluster_manager.h"
#include <algorithm>
#include <cmath>
#include "rmsd.h"

#ifdef WITH_MPI

#include "mpi.h"
#include "mpicoordinator.h"
#include "mpimessaging.h"

#endif

/*
 * We just need the protein collection to initialize an instance
 */
cluster_manager::cluster_manager(vector<vector<atom> >* proteins)
{
	this->population = NULL; // if the proteins are provided directly there will be no population
	this->proteins = proteins;
	
	this->initialize_distances(proteins->size());
}

/*
 * This constructor is used when processing clustering information
 * distributedly. The main rank will hold an instance of cluster_manager
 * that doesn't compute all of the rmsd values. It will only serve as a coordinator
 * of rmsd computations that, once concluded, will trigger updates on the distance
 * vectors stored by this instance
 */
cluster_manager::cluster_manager(vector<md_edge_set>* population)
{
	this->population = population;
	this->proteins = NULL; // no specific atomic structures should be used
	
	initialize_distances(population->size());
}

/*
 * Initialize as many cluster_neighborhoods as complexes we have; each of them 
 * will have the prediction number it corresponds to in order to be able to identify it
 * later on when we sort them (from the one with the most neighbors to the least # of neighbors)
 */
void cluster_manager::initialize_distances(size_t num_complexes)
{
	for(size_t complex_index = 0; complex_index < num_complexes; complex_index++)
	{
		distances.push_back(cluster_neighborhood(complex_index));
		unsorted.push_back(cluster_neighborhood(complex_index));
	}
}

/*
 * Intuitively the clustering procedure computes a matrix of size distances.size() X distances.size()
 * However, because the RMSD is a metric we only need to compute one half of the matrix (upper triangular
 * matrix only). If we label the positions of the upper triangular matrix (ignoring diagonals)
 * we have a total of n(n-1)/2 logical cells that need to be computed
 *
 * The two logical index parameters specify the start and end (inclusive) positions whose
 * distances are stored in the double array provided.
 * This method will compare those distances to the threshold and will update the internal neighbor
 * collections accordingly
 */
void cluster_manager::update_neighbors(size_t logical_start_index, size_t logical_end_index, double* distances, double threshold)
{
	// get the row/column position of the first element
	size_t base_index, compared_to_index;
	size_t distance_index = 0; //points to the current distance value analyzed
	size_t n = this->population->size();
	cluster_manager::translate_logical_to_matrix_position(logical_start_index, this->population->size(),
													&base_index, &compared_to_index);

	// analyze joining elements using the distances
	while(logical_start_index <= logical_end_index) {
		// compare the distance to the threshold and update the neighbors
		if(distances[distance_index] <= threshold)
		{
			this->join(base_index, compared_to_index);
		}

		// update all the actual matrix indices and the logical one
		compared_to_index++;
		if(compared_to_index == n) { // move to the next row
			base_index++;
			compared_to_index = base_index + 1;
		}
		
		logical_start_index++;
		distance_index++;
	}
}

/*
 * Calculate the RMSD values in a region of the distances matrix
 * identified by two logical indices.
 *
 * The distances variable should point to an appropriately allocated array
 * that can hold all distances calculated
 */
void cluster_manager::calculate_distances(size_t logical_start_index, size_t logical_end_index,
											vector<pdb>* pdbs, vector< vector<transformations*> >* predictions,
											double* distances)
{
	vector<atom> base, compared_to;
	size_t base_index, compared_to_index;
	size_t distance_index = 0;
	size_t n = this->population->size();
	cluster_manager::translate_logical_to_matrix_position(logical_start_index, n,
															&base_index, &compared_to_index);
															
	while(logical_start_index <= logical_end_index) {
		base = this->population->at(base_index).get_transformed_complex(*pdbs, *predictions);
		compared_to = this->population->at(compared_to_index).get_transformed_complex(*pdbs, *predictions);
		
		//store the distance
		distances[distance_index++] = calculate_allatom_rmsd(base, compared_to);

		// update all the actual matrix indices and the logical one
		compared_to_index++;
		if(compared_to_index == n) { // move to the next row
			base_index++;
			compared_to_index = base_index + 1;
		}
		
		logical_start_index++;
	}
}

/*
 * Intuitively the clustering procedure computes a matrix of size distances.size() X distances.size()
 * However, because the RMSD is a metric we only need to compute one half of the matrix (upper triangular
 * matrix only). If we label the positions of the upper triangular matrix (ignoring diagonals)
 * we have a total of n(n-1)/2 logical cells that need to be computed
 * This method transforms a logical index into the actual row/column position
 */
void cluster_manager::translate_logical_to_matrix_position(size_t logical_index, size_t n,
													size_t* base_index_out, size_t* compared_to_index_out)
{
	size_t skipped;
	size_t base_index = 0;
	size_t elements_in_row = n - 1; //elements in row 0
	for(skipped = 0; logical_index >= (skipped + elements_in_row); ) {
		base_index++; // go to next row
		skipped += elements_in_row; //we're skipping all elements in this row
		elements_in_row--; // the next one has one less element
	}
	
	*base_index_out = base_index;
	
	// this is the column position
	*compared_to_index_out = (n - elements_in_row) + // + index of first logical element
									(logical_index - skipped); // logical offset
}

/*
 * Inverse of the previous operation, namely, convert from row/column values
 * to the logical view of the rmsd's computed
 */
size_t cluster_manager::translate_matrix_to_logical_position(size_t base_index, size_t compared_to_index, size_t n)
{
	size_t logical_index = 0;
	size_t elements_in_row = n - 1; //elements in row 0
	while(base_index--) {
		logical_index += elements_in_row;
		elements_in_row--;
	}
	
	logical_index += compared_to_index - (n - elements_in_row);
	
	return logical_index;
}

/*
 * Triggers the start of the clustering process, using as "clustering
 * distance" . After this method is called, it will be possible to call
 * the appropriate method to retrieve the clustering results. After this function
 * exits, the internal distance-management structure will be sorted from the
 * complex with the most neighbors (first position) to the one with the least
 * neighbors (last position)
 */
void cluster_manager::cluster(double threshold)
{
	vector<atom> base, compared_to;
	size_t total = proteins->size();
	for(size_t base_index = 0; base_index < total; base_index++)
	{
		base = proteins->at(base_index);

		//for(size_t compared_to_index = 0; compared_to_index < total; compared_to_index++)
		for(size_t compared_to_index = base_index + 1; compared_to_index < total; compared_to_index++)
		{
			/*
			if(base_index == compared_to_index)
			{
				continue; // skip if it's the same one
			}*/
	
			compared_to = proteins->at(compared_to_index);

			double rmsd = calculate_allatom_rmsd(base, compared_to);

			// If it's within the threshold add it as a neighbor of the base prediction
			if(rmsd <= threshold)
			{
				this->join(base_index, compared_to_index);
			}
		}
	}

}

/*
 * Updates the internal neighbor structures in order to
 * joing these two elements as neighbors
 */
void cluster_manager::join(size_t base_index, size_t compared_to_index)
{
	distances[base_index].add_neighbor(compared_to_index);
	unsorted[base_index].add_neighbor(compared_to_index);
	// and we should also add the symmetric relation
	distances[compared_to_index].add_neighbor(base_index);
	unsorted[compared_to_index].add_neighbor(base_index);
}

/*
 * This static clustering method has the sole purpose of coordinating
 * mpi distributed clustering if it's available. Internally, it creates a
 * cluster_manager instance and calls its methods
 */
vector<cluster_neighborhood> cluster_manager::cluster_population(vector<md_edge_set>* population,
																vector<pdb>* pdbs,
																vector< vector<transformations*> >* predictions,
																double threshold)
{

#ifdef WITH_MPI
	int numtasks;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	
	// Calculate sizes for message passing purposes
	int auxranks = numtasks - 1;
	size_t n = population->size();
	size_t total_distances = (n * (n - 1)) / 2;
	size_t distances_per_rank = floor(total_distances / ((float)auxranks));
	// if the number of distance calculations is not a multiple of auxranks
	// we will have some pending calculations that will be performed by the main rank
	size_t remaining_distances = total_distances % auxranks;
	
	//allocate a buffer that can be used to receive the distance updates from the auxiliary ranks
	size_t distances_capacity = remaining_distances > distances_per_rank ? remaining_distances : distances_per_rank;
	double* distances = new double[distances_capacity];

	if(distances_per_rank > 0) {
		// Before clustering we need to send a population update message to the auxiliary ranks
		mpicoordinator::send_distributed_population_update_request(population);
	}
	
	// create a cluster manager that uses a population instead of the explicit atomic structures
	cluster_manager cluster_handler(population);
	
	if(distances_per_rank > 0) {
		// 1. send a message to all auxiliary ranks that they should cluster their expected portion
		// of the distances matrix (each rank knows based on the population size and rank id)
		mpicoordinator::send_distributed_clustering_request(pdbs->size() - 1);
	}
	
	// 2. Before receiving the results check if there are distances left to compute
	// this happens because each rank will compute floor(total / n) and the remainder will be
	// computed by this rank (the main one)
	// In general the remainder could be considerably less than the rest so we manage to use the CPU
	// while the rest of the ranks finish their computations
	if(remaining_distances) {
		size_t logical_start_index = total_distances - remaining_distances; // index of first one pending
		size_t logical_end_index = total_distances - 1; // index of the actual last element
		
		cluster_handler.calculate_distances(logical_start_index, logical_end_index,
												pdbs, predictions, distances);
		// 2.1 Update the neighbors with the distances computed by the main rank
		cluster_handler.update_neighbors(logical_start_index, logical_end_index, distances, threshold);	
	}
	
	// 3. Receive one message per numtask and use the distances they provide to update the neighbors
	while(auxranks-- && distances_per_rank > 0) {
		int error_code = MPI_Recv(distances, distances_per_rank, MPI_DOUBLE, MPI_ANY_SOURCE,
									MPI_MESSAGE_CLUSTER_RECEIVE_TAG, MPI_COMM_WORLD, &status);
		int rank_number = status.MPI_SOURCE;
			
		if(error_code != MPI_SUCCESS) {
			cout << "Error detected while receiving clustering info from auxiliary process" << endl;
		}
		// infer the indices based on the rank that sent the message and update the neighbors
		size_t logical_start_index = distances_per_rank * (rank_number - 1);
		size_t logical_end_index = logical_start_index + (distances_per_rank - 1); //inclusive
		
		cluster_handler.update_neighbors(logical_start_index, logical_end_index, distances, threshold);
	}
	
	delete [] distances;
	// 4. At this point all the neighbors should be correctly set and we can compute the cluster centers
	
#else
	// gather the decoys atoms in order to proceed with the clustering
	vector<vector<atom> > decoys;
	for(size_t decoy_index = 0; decoy_index < population->size(); decoy_index++)
	{
		decoys.push_back(population->at(decoy_index).get_transformed_complex(*pdbs, *predictions));
	}

	/* perform the clustering */
	cluster_manager cluster_handler(&decoys);
	cluster_handler.cluster(threshold);
#endif
	// true means that we will use the highest ranked element of a cluster as its "center"
	vector<cluster_neighborhood> centers = cluster_handler.get_cluster_centers(true);
	return centers;

}

/*
 * Assuming that "cluster" has been called already, this method
 * analyzes the distances calculated and returns a vector with the information
 * corresponding to the cluster centers
 * If the parameter passed is false, the function will consider as centers the elements
 * that have the most number of neighbors. However, if the parameter passed is true
 * the function will first identify the element with the most neighbors too but, it will
 * return whichever neighbor has the lowest (better) rank. This means that although a prediction
 * could have the highest number of neighbors, if one of those neighbors is highly ranked,
 * that will be the one kept
 */
vector<cluster_neighborhood> cluster_manager::get_cluster_centers(bool is_top_ranked_center)
{
	/* The first cluster center will be the complex with the most neighbors and all of it's neighbors
	 * will be skipped in order to find the second most populated cluster (without the neighbors of the
	 * first cluster). The process is repeated until we run out of elements
	 */
	// To do that we need the distances to be sorted
	sort(distances.begin(), distances.end(), cluster_neighborhood::compare);

	// This structure is used to flag the neighbors that must be skipped
	vector<bool> skipped(distances.size(), false);

	// We create a new collection to return, because not all of the elements in distance will
	// be part of the result
	vector<cluster_neighborhood> centers;

	// current_center is the one that we'll add to the result
	cluster_neighborhood current_center;
	// ranking corresponding to the center
	size_t current_center_rank;

	for(size_t distance_index = 0; distance_index < distances.size(); distance_index++)
	{
		current_center = distances[distance_index];
		current_center_rank = current_center.get_prediction_number();
		if(!skipped[current_center_rank]) // if it hasn't been marked as deleted, add the info
		{
			// mark this one in order to skip it if it's found in the future
			skipped[current_center_rank] = true;


			// mark all the neighbors as "to be skipped"
			vector<size_t>* neighbors = distances[distance_index].get_neighbors();
			for(size_t neighbor_index = 0; neighbor_index < neighbors->size(); neighbor_index++)
			{
				// what we store as neighbor is actually the prediction rank
				size_t current_neighbor_rank = neighbors->at(neighbor_index);
				// if we want to keep the top ranked element then compare the current against the ranks of neighbors
				if(is_top_ranked_center && (!skipped[current_neighbor_rank]) // we haven't skipped this neighbor
					&& current_neighbor_rank < current_center_rank) // and the rank of the neighbor is better
				{
					current_center = unsorted[current_neighbor_rank];
					current_center_rank = current_neighbor_rank;
				}
				// regardless of what happens we want to leave the mark
				skipped[current_neighbor_rank] = true;

			}
				
			centers.push_back(current_center);
		}
	}

	return centers;
}

