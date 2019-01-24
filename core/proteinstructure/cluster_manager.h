#ifndef _CLUSTER_MANAGER_H_
#define _CLUSTER_MANAGER_H_

#include <vector>
#include "atom.h"
#include "cluster_neighborhood.h"
#include "mdgraph.h"

using std::vector;

/*
 * Instances of this class can cluster a collection of proteins, each of them
 * represented as a vector<atom>. It will calculate all-against-all pairwise RMSD
 * and, for each protein, will count the number of neighbors, in terms of the RMSD
 * threshold provided
 */
class cluster_manager
{
	private:
		/*
		 * Each vector<atom> represents one protein complex that needs to be clustered
		 */
		vector<vector<atom> >* proteins;
		/*
		 * Alternate source of information for clustering (used when MPI is available)
		 */
		vector<md_edge_set>* population;
		/*
		 * Structure that holds the distances between each protein complex
		 * The size of this vector should be equal to the size of proteins
		 */
		vector<cluster_neighborhood> distances;
		/* This is a copy that maintains the original order before sorting in cluster() */
		vector<cluster_neighborhood> unsorted;

		/*
		 * Initialize as many cluster_neighborhoods as complexes we have; each of them 
		 * will have the prediction number it corresponds to in order to be able to identify it
		 * later on when we sort them (from the one with the most neighbors to the least # of neighbors)
		 */
		void initialize_distances(size_t num_complexes);
		/*
		 * Intuitively the clustering procedure computes a matrix of size distances.size() X distances.size()
		 * However, because the RMSD is a metric we only need to compute one half of the matrix (upper triangular
		 * matrix only). If we label the positions of the upper triangular matrix (ignoring diagonals)
		 * we have a total of n(n-1)/2 logical cells that need to be computed
		 * This method transforms a logical index into the actual row/column position
		 */
		static void translate_logical_to_matrix_position(size_t logical_index, size_t n,
															size_t* base_index_out, size_t* compared_to_index_out);
		/*
		 * Inverse of the previous operation, namely, convert from row/column values
		 * to the logical view of the rmsd's computed
		 */
		static size_t translate_matrix_to_logical_position(size_t base_index, size_t compared_to_index, size_t n);
	public:
		/*
		 * We just need the protein collection to initialize an instance
		 */
		cluster_manager(vector<vector<atom> >* proteins);
		/*
		 * This constructor is used when processing clustering information
		 * distributedly. The main rank will hold an instance of cluster_manager
		 * that doesn't compute all of the rmsd values. It will only serve as a coordinator
		 * of rmsd computations that, once concluded, will trigger updates on the distance
		 * vectors stored by this instance
		 */
		cluster_manager(vector<md_edge_set>* population);
		
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
		void update_neighbors(size_t logical_start_index, size_t logical_end_index, double* distances, double threshold);
		/*
		 * Calculate the RMSD values in a region of the distances matrix
		 * identified by two logical indices.
		 *
		 * The distances variable should point to an appropriately allocated array
		 * that can hold all distances calculated
		 */
		void calculate_distances(size_t logical_start_index, size_t logical_end_index,
													vector<pdb>* pdbs, vector< vector<transformations*> >* predictions,
													double* distances);

		/*
		 * Triggers the start of the clustering process, using as "clustering
		 * distance" . After this method is called, it will be possible to call
		 * the appropriate method to retrieve the clustering results. After this function
		 * exits, the internal distance-management structure will be sorted from the
		 * complex with the most neighbors (first position) to the one with the least
		 * neighbors (last position)
		 */
		void cluster(double threshold);
		/*
		 * Updates the internal neighbor structures in order to
		 * joing these two elements as neighbors
		 */
		void join(size_t base_index, size_t compared_to_index);
		/*
		 * This static clustering method has the sole purpose of coordinating
		 * mpi distributed clustering if it's available. Internally, it creates a
		 * cluster_manager instance and calls its methods
		 * It will return a vector with the clustering results, the same way get_cluster_centers does
		 */
		static vector<cluster_neighborhood> cluster_population(vector<md_edge_set>* population,
																vector<pdb>* pdbs,
																vector< vector<transformations*> >* predictions,
																double threshold);


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
		vector<cluster_neighborhood> get_cluster_centers(bool is_top_ranked_center);
};

#endif
