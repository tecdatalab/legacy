#include "mdgraph.h"
#include "utils.h"

using namespace std;

#include "contact.h"
#include "scoring.h"
#include "mpicoordinator.h"

MTRand rng;

/* Definition of the scoring scheme variable */
score_type_t md_edge_set::scoring_scheme;

/*
 * Generate a graph that will be used as base to create the initial population.
 * The md_graph instance contains a vertex for each value from 0 up to
 * number_of_vertices - 1. It's also going to be a complete graph, that will be used as a base to
 * generate new candidates
 */
md_graph generate_base_graph(int number_of_vertices,
    vector<pair<md_vertex,md_vertex> > interactions)
{
	md_vertex_set vertices;
	// create the number of vertices we need, labeled from 0
	// up to number_of_vertices - 1
	for(int i = 0; i < number_of_vertices; i++)
	{
		vertices.push_back(i);
	}

/// REMOVE
cout << "++++Allowed interactions\n";
for(size_t inter = 0; inter < interactions.size(); inter++) {
  cout << interactions[inter].first << " " <<interactions[inter].second<<endl;
}
cout << "++++END\n";
	
	// create all edges that represent the conformations that were
	// actually generated using ZDOCK or LZerD
	md_edge_set edges;
	for(int i = 0; i < number_of_vertices; i++)
	{
		for(int j = i + 1; j <	number_of_vertices; j++)
		{
      // If interactions are limited, check that the pair is included in them.
      pair<md_vertex,md_vertex> edge_pair = make_pair(i, j);
      if(interactions.empty() || // empty when no restrictions
         std::find(interactions.begin(), interactions.end(), edge_pair) !=
             interactions.end()) {
			  // since this is the base one
        // it is not necessary to set a random prediction number
			  md_edge current_edge(i, j, 0);
			  edges.add_edge(current_edge);
      }
		}
	}
	
	md_graph base_graph(vertices, edges);
	
	return base_graph;
}

/*
 * Creates population_size amount of edge_sets which will be used as the
 * Generation 0 population for the GA process
 */
vector<md_edge_set> generate_initial_population(vector< vector<transformations*> > all_predictions, md_edge_set edges, int population_size)
{
	// set random seed
	//srand ( time(NULL) );
	vector<md_edge_set> population;
	population.reserve(population_size);
	while(population_size--)
	{
		population.push_back(generate_individual(all_predictions, edges));
	}
	return population;
}

/*
 * This function assumes that the stream is already placed on the first line that represents a
 * candidate (after the headers) and creates a md_edge_set instance for each line
 */
vector<md_edge_set> load_population(ifstream& md_file)
{
	vector<md_edge_set> population;
	string edges_str;
	int clashes;
	double pairwise_score, score;
	//each line contains 3 elements: the edges descriptor, the zdock score and the clashes
	// assume we have at least one
	md_file >> edges_str;
	md_file >> pairwise_score;
	md_file >> clashes;
	md_file >> score;

	do
	{
		md_edge_set individual;
		// update the score values for the candidate, specified on the file
		individual.set_pairwise_score(pairwise_score);
		individual.set_score(score);
		individual.set_clashes(clashes);
		vector<string> edge_set_str = split(edges_str, ';');
		// Go through the string representation of each edge and extract the id of the
		// receptor, ligand and the prediction number used
		for(size_t edge_index = 0; edge_index < edge_set_str.size(); edge_index++)
		{
			string edge_str = edge_set_str[edge_index];
			vector<string> edge_info = split(edge_str, ',');
			md_edge new_edge;
			new_edge.receptor = atoi(edge_info[0].c_str());
			new_edge.ligand = atoi(edge_info[1].c_str());
			new_edge.prediction_number = atoi(edge_info[2].c_str());
			individual.add_edge(new_edge);
		}
		population.push_back(individual);
		// update for next iteration
		md_file >> edges_str;
		md_file >> pairwise_score;
		md_file >> clashes;
		md_file >> score;
	} while(!md_file.eof());

	return population;
}

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
              DisjointSets& vertex_sets, bool new_prediction_number)
{
	md_edge_set available_edges = base_edges;
	bool edge_added = false; // we will loop until we can add a valid edge
	while(!edge_added)
	{
		// select a random candidate to add
		int rand_position = (int)floor(rng.randExc(available_edges.size()));
		md_edge candidate = available_edges[rand_position];
		// and remove it because we're not considering it anymore
		available_edges.erase(available_edges.begin() + rand_position);
		// get the set that each vertex belongs to
		int receptor_set = vertex_sets.FindSet(candidate.receptor);
		int ligand_set = vertex_sets.FindSet(candidate.ligand);
		// if they're not connected yet
		if(receptor_set != ligand_set)
		{
			// set the prediction number used if requested
			if(new_prediction_number)
			{
				transformations* prediction_set = all_predictions[candidate.receptor][candidate.ligand];
				int rand_prediction_num = (int)floor(rng.randExc(prediction_set->get_size()));
				candidate.prediction_number = rand_prediction_num;
			}
			// add the edge
			target.add_edge(candidate);
			// join the vertex sets
			vertex_sets.Union(receptor_set, ligand_set);
			edge_added=true;
		}
		// else: repeat the loop, continue looking
	}
}

/*
 * Create a completely new individual based on the predictions and edges provided
 */
md_edge_set generate_individual(vector< vector<transformations*> > all_predictions, md_edge_set base_edges)
{
	md_edge_set individual;
	// the number of vertices is equal to the number of columns in the predictions matrix
	int total_vertices = all_predictions.size();
	DisjointSets vertex_set(total_vertices);
	// NULL indicates that we don't have a previously create DisjointSets instance to help
	// true indicates that prediction numbers for each edge will be selected randomly
	complete_graph(individual, all_predictions, base_edges, vertex_set, true);
	return individual;
}

/*
 * Given an existing indidividual, a sert of predictions and source edges,
 * add as many edges to the inidividual as necessary in order to create a
 * spanning tree.
 */
void complete_graph(md_edge_set& individual, vector< vector<transformations*> >& all_predictions, md_edge_set& base_edges,
						DisjointSets& s, bool new_prediction_number)
{
	// copy the base edges because they are going to be successively removed one by one
	// until we get a complete spanning tree
	md_edge_set available_edges = base_edges;
	// the number of vertices is equal to the number of columns in the predictions matrix
	size_t total_vertices = all_predictions.size();
	while(individual.size() < (total_vertices - 1))
	{
		// add random edges until we connect all of them
		add_edge(all_predictions, available_edges, individual,
		         	s, new_prediction_number);
	}
}


/*
 * Create a new edge_set by removing one of the edges from the original individual
 * and later reconnecting the graph, by randomly selecting edges from base
 */
md_edge_set mutate(md_edge_set individual, md_edge_set& base, vector< vector<transformations*> >& all_predictions)
{
	md_edge_set mutated_edge_set;
	// select one of the edges to be deleted
	size_t deleted_edge = (size_t)floor(rng.randExc(individual.size()));
	// the number of vertices is equal to the number of columns in the predictions matrix
	int total_vertices = all_predictions.size();
	// Since we will delete an edge and try to reconnect the graph, we need
	// a DisjointSets data structure
	DisjointSets vertex_sets(total_vertices);
	for(size_t edge_index = 0; edge_index < individual.size(); edge_index++)
	{
		if(edge_index == deleted_edge)
		{
			// don't add the edge we're discarding
			continue;
		}
		md_edge current_edge = individual[edge_index];
		mutated_edge_set.add_edge(current_edge);
		
		// get the set that each vertex belongs to
		int receptor_set = vertex_sets.FindSet(current_edge.receptor);
		int ligand_set = vertex_sets.FindSet(current_edge.ligand);
		// if they're not connected yet in the disjoint sets structures
		if(receptor_set != ligand_set)
		{
			// join the vertex sets
			vertex_sets.Union(receptor_set, ligand_set);
		}
	}
	// finally, since we deleted one edge, we neeed to re-connect the graph
	// the "true" parameter indicates that we want prediction numbers to be selected
	// randomly (as opposed to keeping the prediction number that the base edges have)
	complete_graph(mutated_edge_set, all_predictions, base, vertex_sets, true);

	return mutated_edge_set;
}

/*
 * Given two sets of edges, create an offspring out of them, by randomly
 * selecting edges from the ones that the parents have, until a connected
 * graph is created.
 */
md_edge_set recombine(md_edge_set parent1, md_edge_set parent2, vector< vector<transformations*> >& all_predictions)
{
	// add parent1's edges...
	md_edge_set parent_combination = parent1;
	// and then copy the rest
	for(size_t parent2_index = 0; parent2_index < parent2.size(); parent2_index++)
	{
		md_edge current_edge = parent2[parent2_index];
		parent_combination.add_edge(current_edge);
	}
	// Create an empty edge set
	md_edge_set child;
	// the number of vertices is equal to the number of columns in the predictions matrix
	int total_vertices = all_predictions.size();
	// Since we will delete and edge and try to reconnect the graph, we need
	// a DisjointSets data structure
	DisjointSets vertex_sets(total_vertices);
	// and call the function that will convert it into a connected graph
	// NULL indicates that we don't have a previously create DisjointSets instance to help
	// false indicates that prediction numbers for each base edge will be kept as is
	complete_graph(child, all_predictions, parent_combination, vertex_sets, false);
	return child;
}

	
/*
 * This operator is overriden because we need to sort edge sets
 * to prune the least fit subjects from each generation.
 * The sort function implemented in the algorithm library, requires this.
 * The function returns true if the first instance has higher zdock_score,
 * or if it has higher clash score, in case of a tie
 */
bool md_edge_set::compare(md_edge_set edge_set1, md_edge_set edge_set2)
{
	// Create two cases depending on the type of scoring used
	if(scoring_scheme == score_type_pairwise)
	{
		// In this case, both zdock and lzerd use a "higher is better" scheme
		// we need to return true if the first element must go first, namely, if it's higher
		// this instance is less than compare to if the zdock score is lower...
		return (edge_set1.pairwise_score > edge_set2.pairwise_score) ||
			// or if they have equal score, but the clashes are less (used as a "tie-breaker")
			(edge_set1.pairwise_score == edge_set2.pairwise_score &&
			 edge_set1.clashes > edge_set2.clashes);
	}
	else // physics
	{
		// In this case, the physics based scheme is "lower is better" so we must return true
		// if the first one must go first, in other words, if the first one is lower
		return (edge_set1.score < edge_set2.score) ||
			// or if they have equal score, but the clashes are less (used as a "tie-breaker")
			(edge_set1.score == edge_set2.score &&
			 edge_set1.clashes > edge_set2.clashes);
	}
}

/*
 * It pushes the new edge into the collection of edges that this instances
 * as well as setting the score update flag to true (so that it is indeed updated
 * the next time that update_score is called)
 */
void md_edge_set::add_edge(md_edge new_edge)
{
	push_back(new_edge);
	needs_score_update = 1;
}

/*
 * Assuming the new values that determine the score of this instance
 * have been already determined, this method just takes the values in the
 * parameter and updates the instance variables accordingly
 */
void md_edge_set::update_score(double new_score, double new_pairwise_score, int new_clashes, soroban_t new_soro_score)
{
	this->score = new_score;
	this->pairwise_score = new_pairwise_score;
	this->clashes = new_clashes;
	this->soro_score.set_base(new_soro_score);
	/* once we finish set the score update flag to false */
	needs_score_update = 0;
}

/*
 * Overload that receives a struct instead of the separate values
 */
void md_edge_set::update_score(md_edge_set_score_t new_score)
{
	this->update_score(new_score.score, new_score.pairwise_score, new_score.clashes, new_score.soro_score.get_base());
}
		
/*
 * Instance version of the methods that uses "this" as the set of edges
 */
void md_edge_set::update_score(vector<pdb>& proteins, vector< vector<transformations*> >& predictions, int clash_cutoff,
	double* weights)
{
	/* since this is a costly operation, only recompute it if there has been a change */
	if(needs_score_update) {
		md_edge_set_score_t new_score = update_score(this, proteins, predictions, clash_cutoff, weights);
		this->update_score(new_score);
	}
}


/*
 * Since adding up the zdock scores and calculating the clash score
 * is something that should be computed as few times as possible, it is
 * necessary that this function is called for all edge_set instances
 * before ranking the elements that belong to a certain generation
 * It is necessary that a clash threshold is provided for efficiency purposes (if it's
 * reached, then the clash calculation will stop). clash_cutoff = 0 means that it should be
 * ignored
 */
md_edge_set_score_t md_edge_set::update_score(vector<md_edge>* all_edges, vector<pdb>& proteins,
		vector< vector<transformations*> >& predictions, int clash_cutoff,
		double* weights)
{
	md_edge_set_score_t result;

	// reset the values
	result.pairwise_score = 0.0;
	result.clashes = 0;
	result.score = 0;
	// calculate pairwise score
	for(size_t i = 0 ; i < all_edges->size(); i++)
	{
		md_edge edge = (*all_edges)[i];
		result.pairwise_score += predictions[edge.receptor][edge.ligand]->get_score(edge.prediction_number);
	}
	// Normalize by the number of edges
	result.pairwise_score /= all_edges->size();
			
	// get the transformed atoms
	vector< vector<atom> > transformed_proteins = get_transformed_atoms(proteins, predictions, all_edges);
	// calculate clashes score
	vector<atom> cumulative_atoms; // start with the first protein, which is the base
	for(size_t protein_index = 0; protein_index < transformed_proteins.size(); protein_index++)
	{
		vector<atom> next_protein = transformed_proteins[protein_index];
		// get the clashes between the current set and the next protein
		update_clash_count(cumulative_atoms, next_protein, result.clashes, clash_cutoff);
		// add the current protein to the current list of atoms
		cumulative_atoms.reserve(cumulative_atoms.size() + next_protein.size());
		cumulative_atoms.insert(cumulative_atoms.begin() + cumulative_atoms.size(),
						next_protein.begin(), next_protein.end());
	}
			
	// update the energy score only if the clash cutoff needs to be ignored (clash_cutoff = 0) or if the clashes are
	// less than the cutoff
	if(clash_cutoff == 0 || result.clashes < clash_cutoff)
	{
		result.soro_score = compute_energy(transformed_proteins, weights);
		result.score = result.soro_score.calculate_weighted_score();
	}

	return result;
}

/*
 * Takes a collection of md_edge_set instances and call their respective
 * update_score methods. Optionally, it can be indicated that only the first
 * n elements will have their scores updated
 */
void md_edge_set::update_score(vector<md_edge_set>* population, vector<pdb>& proteins,
                                vector< vector<transformations*> >& predictions, int clash_cutoff,
				double* weights,int up_to_index)
{
	int n = up_to_index == -1 ? population->size() : up_to_index;
#ifdef WITH_MPI
	MPI_Status status;
	// the messages received will deposit the result here
	mpi_message_update_score_receive_t result_msg;
	// get the number of processes
	int numtasks;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	
	MPI_Datatype main_to_aux_type = mpi_type_update_score_send();
	MPI_Datatype aux_to_main_type = mpi_type_update_score_receive();
	int num_edges = proteins.size() - 1; // this is the total edges in the spanning tree
	mpi_message_update_score_send_t* edges_message = new mpi_message_update_score_send_t[num_edges];

	// These 2 values track the elements sent to the distributed processes and
	// if their respective results were received
	int sent_index = 0;
	int sent_items = 0; // the difference between sent_index and this one is that this is a counter of the messages
				// not the actual index sent

	// send the first messages
	// if there are more processes than elements left we'll stop when we run out of elements
	for(; (sent_index < n) && sent_items < (numtasks - 1); sent_index++) {
		md_edge_set current_edges = (*population)[sent_index];
		if(current_edges.needs_update()) {
			mpicoordinator::send_distributed_score_request((*population)[sent_index], edges_message,
	                                                        main_to_aux_type, sent_index, ++sent_items);
		}
	}
	

	// check from ranks 1 to numtasks, receive the score and send a new message
	// if we still have elements to score
	if(sent_items > 0) {
		do {
			// make a copy of the expected number of messages back
//			int currently_pending = sent_items;
			// and reset the sent counter
//			sent_items = 0;
//			for(int rank_number = 1; rank_number <= currently_pending; rank_number++) {
//				MPI_Recv(&result_msg, 1, aux_to_main_type, rank_number, MPI_MESSAGE_UPDATE_SCORE_RECEIVE_TAG, MPI_COMM_WORLD, &status);
				int error_code = MPI_Recv(&result_msg, 1, aux_to_main_type, MPI_ANY_SOURCE,
								MPI_MESSAGE_UPDATE_SCORE_RECEIVE_TAG, MPI_COMM_WORLD, &status);
				int rank_number = status.MPI_SOURCE;
			
				if(error_code != MPI_SUCCESS) {
					cerr << "Error detected while receiving update score message on auxiliary process" << endl;
				}
				// update the number of items that are pending
				sent_items--;
				
				// and update the score
				
				(*population)[result_msg.population_index].update_score(result_msg.score, result_msg.pairwise_score,
											result_msg.clashes, result_msg.soro_score);
			
				// at this point sent_index refers to the last item sent, we move it to the next one
				// that needs sending (because it needs update
				// if we have pending tasks, send a new message to this rank
				while(sent_index < n && (!(*population)[sent_index].needs_update())) {
					sent_index++;
				}

				// if, after skipping over, we still have something to send then do so
				if(sent_index < n) {
					// increment the sent_items counter using prefix notation, because that
					// way it will match the rank number
					// sent_index is incremented later because it's already pointing to
					// the element we want to send now, but we want to advance that for the next iteration
//					md_edge_set::send_distributed_score_request((*population)[sent_index], edges_message,
//	                                			                        main_to_aux_type, sent_index, ++sent_items);
					mpicoordinator::send_distributed_score_request((*population)[sent_index], edges_message,
	                                			                        main_to_aux_type, sent_index, rank_number);
					sent_items++; // one more pending
			                sent_index++;
				}
//			}
		} while(sent_items);
	}

	MPI_Type_free(&main_to_aux_type);
	MPI_Type_free(&aux_to_main_type);
	delete [] edges_message;
#else
	for(int population_index = 0; population_index < n; population_index++) {
		(*population)[population_index].update_score(proteins, predictions, clash_cutoff, weights);
	}
#endif
}

/*
 * Two edge sets are considered equal if they have the same size and
 * they have exactly the same edges. This method returns true if this
 * instance is equal to the one provided as parameter
 */
bool md_edge_set::equals(md_edge_set other_set)
{
	bool matches;
	if(other_set.size() != this->size())
	{
		return false;
	}
	else
	{
		for(size_t edge_index = 0; edge_index < this->size(); edge_index++)
		{
			md_edge current = (*this)[edge_index];
			md_edge matched_edge = other_set.find_edge(current.receptor, current.ligand, &matches);
			/* if just one of them is different return false 
			 * also, find_edge doesn't check the prediction number so, if we found a receptor/ligand
			 * match we need to check if the prediction number is different
			 */
			if(!matches || (matches && matched_edge.prediction_number != current.prediction_number))
			{
				return false;
			}
		}
		/* if we get to the end all of them matched */
		return true;
	}
}

/*
 * Peform a linear search for the edge that contains the receptor/ligand
 * connection specified.
 */
md_edge md_edge_set::find_edge(vector<md_edge>* all_edges, md_vertex rec, md_vertex lig, bool* found)
{
	md_edge result;
	/* by default we assume it wasn't found */
	(*found) = false;
	for(size_t edge_index = 0; edge_index < all_edges->size(); edge_index++)
	{
		md_edge current = (*all_edges)[edge_index];
		if(current.receptor == rec && current.ligand == lig)
		{
			result = current;
			(*found) = true;
			break;
		}
	}
	return result;
}

/*
 * Invokes get_transformed_atoms and integrates all chains into one single vector with
 * all the atoms after applying the appropriate transformations.
 * After the first time this method is invoked (for this instance), the function will store the result
 * in transformed_complex. Subsequent calls to this will not recalculate the transformation, they will just
 * return the stored vector
 */
vector<atom> md_edge_set::get_transformed_complex(vector<pdb>& proteins, vector< vector<transformations*> >& predictions)
{
	if(transformed_complex.empty()) //calculate the transformed positions
	{
		vector<vector<atom> > all_chains = get_transformed_atoms(proteins, predictions);
		for(size_t protein_index = 0; protein_index < all_chains.size(); protein_index++)
		{
			vector<atom> next_protein = all_chains[protein_index];
			// add the current protein to the current list of atoms
			transformed_complex.reserve(transformed_complex.size() + next_protein.size());
			transformed_complex.insert(transformed_complex.begin() + transformed_complex.size(),
						next_protein.begin(), next_protein.end());
		}
	} // else, if it's already full then the transformation has been calculated already (we don't do anything else)
	return transformed_complex;
}

		
/*
 * Before evaluating the score and calculate the clashes, we need to transform the coordinates
 * of the proteins (except the first one, which doesn't need it).
 * This function outputs all the sets of atoms in their final representation, including the
 * first protein, which will be unaltered. Each position of the returned vector will match the
 * vertex index used by each vertex to identify each protein.
 */
vector< vector<atom> > md_edge_set::get_transformed_atoms(vector<pdb>& proteins,
			                        vector< vector<transformations*> >& predictions,
						vector<md_edge>* all_edges)
{
	// Allocate the vector that will hold the atom sets
	vector< vector<atom> > result(proteins.size());
	vector<atom> base_protein = proteins[0].atoms;
	// add the first protein unaltered in position 0
	result[0] = base_protein;
	// for all the rest, we should have a series of rotations that need to be applied
	vector< vector<md_vertex> > rotations = get_rotational_paths(all_edges);
	for(size_t rotation_index=0; rotation_index < rotations.size(); rotation_index++)
	{
		vector<md_vertex> rotation_path = rotations[rotation_index];
		md_vertex protein_index = rotation_path.back();
		vector<atom> cumulative_transformation = proteins[protein_index].atoms;
		vector<atom> after_current_transformation;
		while(rotation_path.size() > 1) // each pair of vertices in the list represents a transformation
		{
			bool found;
			// clear the temporary atom vector used for this iteration's transformation
			after_current_transformation.clear();
			md_vertex lig = rotation_path[rotation_path.size() - 1];
			md_vertex rec = rotation_path[rotation_path.size() - 2];

			md_edge edge = find_edge(all_edges, rec, lig, &found);
			/* it is possible that the edge has the opposite direction */
			if(!found)
			{
				edge = find_edge(all_edges, lig,rec,&found);
				/* Note, this last found value is not checked because one of the two must exist */
			}

			/* the transformation that we need corresponds to the current path segment */
			transformations* current_predictions = predictions[rec][lig];
			current_predictions->transform_atoms(cumulative_transformation, after_current_transformation,
								edge.prediction_number);
			cumulative_transformation = after_current_transformation;
			// remove the last vertex in the path (because we proceed from last to first)
			rotation_path.pop_back();
		}
		result[protein_index] = cumulative_transformation;
	}
	return result;
}

/*
 * Once the edge_set is fixed, we will assume that vertex 0 is the "reference" frame
 * and that we need to rotate all the rest of the proteins to be consistent with protein 0
 * This function returns a list of edges, each representing a transformation, for each of the proteins
 * The result will have N-1 sets of vertices (because protein 0 doesn't need to be moved)
 */
vector< vector<md_vertex> > md_edge_set::get_rotational_paths(vector<md_edge>* all_edges)
{
	// Create a copy of the edges to track which have been visited
	vector<md_edge> these_edges = (*all_edges);
	// The result will be all the rotations that need to be performed for each
	// protein. For example, if we have (0, 2, 4), it means that to obtain
	// the correct conformation of protein 4, we need to rotate it through the
	// edge from 2 to 4 (which represents a transformation) and then the one
	// from 0 to 2
	vector< vector<md_vertex> > result;
	// We assume that vertex 0 is the rotational base
	vector< vector<md_vertex> > candidates_queue;
	vector<md_vertex> initial_base;
	initial_base.push_back(0);//the first candidate path is the one that starts at 0
	candidates_queue.push_back(initial_base);
	do
	{
		// extract one of the candidate paths from the queue
		vector<md_vertex> current_base = candidates_queue.back();
		// and delete it because we won't consider it anymore
		candidates_queue.pop_back();
		// Look for the edges that start with the last element of the current base
		md_vertex current_vertex = current_base.back();
		// 
		// we go through the elements backwards because we could delete elements
		for(int i=these_edges.size()-1; i >= 0; i--)
		{
			md_edge current_edge = these_edges[i];
			/* we need to look for any edges that either start or end with current_vertex,
			 * because those are the ones that will expand the path */
			if(current_edge.receptor == current_vertex || current_edge.ligand == current_vertex)
			{
				// add this path to the result
				md_vertex other_edge = current_edge.receptor == current_vertex ?
							current_edge.ligand : current_edge.receptor;
				vector<md_vertex> correct_path = current_base;//this makes a copy
				correct_path.push_back(other_edge);
				result.push_back(correct_path);

				// but also add 2 paths to the queue, because other vertices could
				// be connected either to the path that initiated this one (branching in other direction)
				// or to correct_path too, in case that there are other paths that don't branch, but rather
				// extend this one
				candidates_queue.push_back(correct_path);
				candidates_queue.push_back(current_base);
				// remove the edge from these_edges because it has been visited
				these_edges.erase(these_edges.begin() + i);
			}
		}
	} while(!candidates_queue.empty());
	return result;
}

/*
 * Sends a string representation of this object to standard output
 */
void md_edge_set::print(ostream& output_stream)
{
	// Format for each decoy: rec_index,lig_index,pred_num rec_index,lig_index,pred_num etc
	// each triple is separated by a space character
	for(size_t i = 0 ; i < this->size(); i++)
	{
		md_edge c = (*this)[i];
		if(i != 0)
		{
			// print a separator because there was an edge before this one
			output_stream << ";";
		}
		output_stream  << c.receptor << "," << c.ligand <<
			"," << c.prediction_number;
	}
	// Add the score and the clashes separated by tabs
	output_stream << "\t" << pairwise_score << "\t" << clashes << "\t" << score << endl;
}

/*
 * Writes this instance to a PDB file
 */
void md_edge_set::write(vector<pdb>& proteins,vector< vector<transformations*> >& predictions,
           string fname, int k)
{
	vector< vector<atom> > transformed_proteins = get_transformed_atoms(proteins, predictions);
	write_complex(fname, k, transformed_proteins);
}
