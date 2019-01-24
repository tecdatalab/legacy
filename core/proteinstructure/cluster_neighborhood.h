#ifndef _CLUSTER_NEIGHBORHOOD_H_
#define _CLUSTER_NEIGHBORHOOD_H_

using namespace std;

/*
 * Small data structure used by the clustering program
 * to store the neighbors (in terms of distance) of each prediction
 */
class cluster_neighborhood
{
	private:
		// Indicates that this cluster_neighborhood is associated to this prediction number
		// (i.e. the prediction with rank N in the original LZerD file)
		size_t this_prediction;
		// Rank of the predictions that are close to this one (RMSD < cutoff)
		vector<size_t> prediction_indices;
	public:
		/*
		 * Default constructor that should only be used to declare variables
		 * that will be later overwritten by assigning some other cluster_neighborhood's
		 * values
		 */
		cluster_neighborhood()
		{
		}
		/*
		 * To create a new neighborhood we need to specify the prediction that
		 * it corresponds to
		 */
		cluster_neighborhood(size_t prediction_number)
		{
			this_prediction = prediction_number;
		}
		/*
		 * Inserts one more neighbor
		 */
		void add_neighbor(size_t prediction_index)
		{
			prediction_indices.push_back(prediction_index);
		}

		/*
		 * Returns the rank of the prediction that this neighborhood represents
		 */
		size_t get_prediction_number()
		{
			return this_prediction;
		}

		/*
		 * Returns the number of neighbors that have been added to this instance
		 */
		size_t get_neighbor_count()
		{
			return prediction_indices.size();
		}

		/*
		 * Returns a collection with the rank of the predictions that are in
		 * this neighborhood
		 */
		vector<size_t>* get_neighbors()
		{
			return &prediction_indices;
		}

		/*
		 * This operator is overriden because we need to sort the neighborhoods
		 * by the number of close predictions to them.
		 * The sort function implemented in the algorithm library, requires this.
		 */
		static bool compare(cluster_neighborhood left, cluster_neighborhood right)
		{
			// We need to return tru if the first element goes first (i.e. if it has more neighbors)
			return left.get_neighbor_count() > right.get_neighbor_count();
		}
};

#endif
