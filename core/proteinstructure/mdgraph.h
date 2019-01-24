#ifndef _MDGRAPH_H_
#define _MDGRAPH_H_

#include "transformations.h"
#include "DisjointSets.h"
#include "MersenneTwister.h"
#include "pdb.h"
#include "soroban_score.h"
#include "md.h"
#include <time.h>
#include <vector>
#include <fstream>

using std::string;
using std::vector;

#define CLASH_ADJUSTMENT 10.0
#define MUTATION_PROBABILITY 0.1
#define RECOMBINATION_PROBABILITY 0.5

#ifdef WITH_MPI
#include "mpi.h"
#include "mpimessaging.h"
#endif

// Each protein chain will be identified by an index. For example, if we had
// chains A, B, C we could assign index 0 to A, index 1 to B and 2 to C
typedef int md_vertex;

// Simply a collection of int's that represent the vertices
typedef vector<md_vertex> md_vertex_set;

/*
 * Mirrors the instance values held by each md_edge_set
 * instance. This was added because of MPI message passing
 */
typedef struct md_edge_set_score_t
{
	double score;
	double pairwise_score;
	int clashes;
	soroban_score soro_score;
} md_edge_set_score_t;

// An edge will be represented as a pair of integers. Each one will
// be an index to the vertex that each end-point of the edge identifies.
// Additionally, they will be characterized as receptor and ligand,
// in order to know which one is fixed and which one was rotated.
// Finally, each edge will contain a prediction number that will determine
// which of the several decoys between receptor and ligand will be used
class md_edge
{
	public:
		md_vertex receptor, ligand;
		int prediction_number;

		md_edge()
		{
		}
	
		md_edge(md_vertex rec, md_vertex lig, int pred_num)
		{
			receptor = rec;
			ligand = lig;
			prediction_number = pred_num;
		}

		md_edge(const md_edge &edge)
		{
			receptor = edge.receptor;
			ligand = edge.ligand;
			prediction_number = edge.prediction_number;
		}
};


// A collection of vertex associations
class md_edge_set: public vector<md_edge>
{
	private:
		double score;
		double pairwise_score;
		int clashes;
		// this variable is false by default and is set to true after an update score
		// is performed on this instance. Each time a new edge is added to the instance
		// it will be set to false again
		bool needs_score_update;
		soroban_score soro_score;
		/*
		 * stores all the atoms (from all chains) after applying the transformation using the information in this graph
		 * This is basically the same than the result returned by the get_transformed_atoms method but, instead of
		 * returning one vector<atom> for each chain, it integrates all the atoms into one vector
		 */
		vector<atom> transformed_complex;
	public:
		// This static variable is set by the main program in order to
		// manage the sorting type used to create the ranking in each generation
		static score_type_t scoring_scheme;
		
		// default constructor, sets everything to zero
		md_edge_set()
		{
			score = pairwise_score = 0;
			clashes = 0;
			needs_score_update = 1;
		}

		/*
		 * Exposes the value of the score update flag
		 */
		bool needs_update()
		{
			return needs_score_update;
		}

		/*
		 * Method used to set the scoring scheme used to sort
		 * the ranking in each generation
		 */
		static void set_scoring_scheme(score_type_t score_type)
		{
			md_edge_set::scoring_scheme = score_type;
		}

		/*
		 * This operator is overriden because we need to sort edge sets
		 * to prune the least fit subjects from each generation.
		 * The sort function implemented in the algorithm library, requires this.
		 * The function returns true if the first instance has higher zdock_score,
		 * or if it has higher clash score, in case of a tie
		 */
		static bool compare(md_edge_set edge_set1, md_edge_set edge_set2);

		/*
		 * It pushes the new edge into the collection of edges that this instances
		 * as well as setting the score update flag to true (so that it is indeed updated
		 * the next time that update_score is called)
		 */
		void add_edge(md_edge new_edge);

		/*
		 * Assuming the new values that determine the score of this instance
		 * have been already determined, this method just takes the values in the
		 * parameter and updates the instance variables accordingly
		 */
		void update_score(double new_score, double new_pairwise_score, int new_clashes, soroban_t new_soro_score);
		/*
		 * Overload that receives a struct instead of the separate values
		 */
		void update_score(md_edge_set_score_t new_score);
		/*
		 * Instance version of the methods that uses "this" as the set of edges
		 */
		void update_score(vector<pdb>& proteins, vector< vector<transformations*> >& predictions, int clash_cutoff,
			double* weights);
	
		/*
		 * Since adding up the zdock scores and calculating the clash score
		 * is something that should be computed as few times as possible, it is
		 * necessary that this function is called for all edge_set instances
		 * before ranking the elements that belong to a certain generation
		 * It is necessary that a clash threshold is provided for efficiency purposes (if it's
		 * reached, then the clash calculation will stop). clash_cutoff = 0 means that it should be
		 * ignored
		 */
		static md_edge_set_score_t update_score(vector<md_edge>* all_edges, vector<pdb>& proteins,
				vector< vector<transformations*> >& predictions, int clash_cutoff,
				double* weights);
	
		/*
		 * Takes a collection of md_edge_set instances and call their respective
		 * update_score methods. Optionally, it can be indicated that only the first
		 * n elements will have their scores updated
		 */
		static void update_score(vector<md_edge_set>* population, vector<pdb>& proteins,
		                                vector< vector<transformations*> >& predictions, int clash_cutoff,
						double* weights,int up_to_index=-1);
	
		/*
		 * Two edge sets are considered equal if they have the same size and
		 * they have exactly the same edges. This method returns true if this
		 * instance is equal to the one provided as parameter
		 */
		bool equals(md_edge_set other_set);
	
		/*
		 * Instance version of find edge that uses "this" as the vector of edges
		 */
		md_edge find_edge(md_vertex rec, md_vertex lig, bool* found) {
			return find_edge(this, rec, lig, found);
		}
		
		/*
		 * Peform a linear search for the edge that contains the receptor/ligand
		 * connection specified.
		 */
		static md_edge find_edge(vector<md_edge>* all_edges, md_vertex rec, md_vertex lig, bool* found);
	
		/*
		 * Invokes get_transformed_atoms and integrates all chains into one single vector with
		 * all the atoms after applying the appropriate transformations.
		 * After the first time this method is invoked (for this instance), the function will store the result
		 * in transformed_complex. Subsequent calls to this will not recalculate the transformation, they will just
		 * return the stored vector
		 */
		vector<atom> get_transformed_complex(vector<pdb>& proteins, vector< vector<transformations*> >& predictions);
	
		/*
		 * Instance version of the method that takes "this" as the edge set
		 */
		vector< vector<atom> > get_transformed_atoms(vector<pdb>& proteins,
					                        vector< vector<transformations*> >& predictions)
		{
			return get_transformed_atoms(proteins, predictions, this);
		}

		/*
		 * Before evaluating the score and calculate the clashes, we need to transform the coordinates
		 * of the proteins (except the first one, which doesn't need it).
		 * This function outputs all the sets of atoms in their final representation, including the
		 * first protein, which will be unaltered. Each position of the returned vector will match the
		 * vertex index used by each vertex to identify each protein.
		 */
		static vector< vector<atom> > get_transformed_atoms(vector<pdb>& proteins,
					                        vector< vector<transformations*> >& predictions,
								vector<md_edge>* all_edges);

		/*
		 * Once the edge_set is fixed, we will assume that vertex 0 is the "reference" frame
		 * and that we need to rotate all the rest of the proteins to be consistent with protein 0
		 * This function returns a list of edges, each representing a transformation, for each of the proteins
		 * The result will have N-1 sets of vertices (because protein 0 doesn't need to be moved)
		 */
		static vector< vector<md_vertex> > get_rotational_paths(vector<md_edge>* all_edges);

		/*
		 * Allows setting the clashes when loading from a file
		 */
		void set_clashes(int new_clashes)
		{
			clashes = new_clashes;
		}

		/*
		 * Allows setting the pairwise score when loading from a file
		 */
		void set_pairwise_score(double new_pairwise_score)
		{
			pairwise_score = new_pairwise_score;
		}

		/*
		 * Allows setting the score when loading from a file
		 */
		void set_score(double new_score)
		{
			score = new_score;
		}

		// Access method that allows read-only access to the score value
		double get_score()
		{
			return score;
		}

		// Access method that allows read-only access to the pairwise score value
		double get_pairwise_score()
		{
			return pairwise_score;
		}

		// Access method that allows read-only access to the clashes value
		double get_clashes()
		{
			return clashes;
		}

		/*
		 * Sends a string representation of this object to standard output
		 */
		void print(ostream& output_stream);
	
		/*
		 * Sends a comma separated list of the soroban score calculated for this complex to
		 * the output stream provided
		 */
		void print_detailed_score(ostream& output_stream)
		{
			soro_score.print(output_stream);
		}

		/*
		 * Writes this instance to a PDB file
		 */
		void write(vector<pdb>& proteins,vector< vector<transformations*> >& predictions,
		           string fname, int k);
};

// A graph will be represented as usual: a collection
// of vertices and edges
class md_graph
{
	public:
		md_vertex_set V;
		md_edge_set E;

		md_graph()
		{
		}
	
		md_graph(md_vertex_set vertices, md_edge_set edges)
		{
			V = vertices;
			E = edges;
		}
};

/*
 * Generate a graph that will be used as base to create the initial population.
 * The md_graph instance returns contains a vertex for each value from 0 up to
 * number_of_vertices - 1, and as many edges as predictions were generated using
 * ZDock
 */
md_graph generate_base_graph(int number_of_vertices,
    vector<pair<md_vertex,md_vertex> > interactions);

/*
 * Creates population_size amount of edge_sets which will be used as the
 * Generation 0 population for the GA process
 */
vector<md_edge_set> generate_initial_population(vector< vector<transformations*> > all_predictions, md_edge_set edges,
													int population_size);

/*
 * Randomly select a new edge (from available_edges) in order to complete the target "graph"
 * (represented by edges in this case). Since in order to capture cycles and graph completeness
 * we use the DisjointSets data structure, it is required that the caller supplies one such instance
 * that contains the appropriate representation of the target graph.
 * The new_prediction_number parameter should be set to true if, once a random edge was selected,
 * it is desired that the prediction number is overwritten with a new one going from 0 to N-1, where
 * N is the number of predictions.
 */
void add_edge(vector< vector<transformations*> > all_predictions, md_edge_set base_edges, md_edge_set& target,
              DisjointSets& vertex_sets, bool new_prediction_number);

/*
 * Create a completely new individual based on the predictions and edges provided
 */
md_edge_set generate_individual(vector< vector<transformations*> > all_predictions, md_edge_set base_edges);

/*
 * Given an existing indidividual, a sert of predictions and source edges,
 * add as many edges to the inidividual as necessary in order to create a
 * spanning tree.
 */
void complete_graph(md_edge_set& individual, vector< vector<transformations*> >& all_predictions, md_edge_set& base_edges,
			DisjointSets& vertex_sets, bool new_prediction_number);

/*
 * Create a new edge_set by removing one of the edges from the original individual
 * and later reconnecting the graph, by randomly selecting edges from base
 */
md_edge_set mutate(md_edge_set individual, md_edge_set& base,
                   	vector< vector<transformations*> >& all_predictions);

/*
 * Given two sets of edges, create an offspring out of them, by randomly
 * selecting edges from the ones that the parents have, until a connected
 * graph is created.
 */
md_edge_set recombine(md_edge_set parent1, md_edge_set parent2,
                      vector< vector<transformations*> >& all_predictions);

/*
 * Get a floating point score for the individual (an edge set) analyzed.
 * It takes into account ZDOCK's score and clashes
 */
double score_individual(md_edge_set& individual, vector<pdb>& proteins, vector< vector<transformations*> >& predictions);
/*
 * This function assumes that the stream is already placed on the first line that represents a
 * candidate (after the headers) and creates a md_edge_set instance for each line
 */
vector<md_edge_set> load_population(ifstream& md_file);

#endif
