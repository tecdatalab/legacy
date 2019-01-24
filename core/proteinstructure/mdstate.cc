#include "mdstate.h"
#include <algorithm>

/*
 * Initializes the simulation state based on the configuration provided
 * Note that this can be either the start of a new simulation or we can
 * be resuming a previous one
 */
mdstate::mdstate(mdga_options* options)
{
	/* first check if we are resuming a previous execution */
	if(options->resume_required())
	{
		read(options->get_resumed_filename());
	}
	else // some of the instance variables are initialized differently if we are not resuming
	{
		this->complex_name = options->get_main_complex_name();
		this->chain_ids = options->get_chain_ids();
		this->current_generation = 0; // 0 because it's the initial population
	}

	/* we don't store the input directory as instance variable but we need it to load the pdb's and predictions */
	string directory = options->get_input_directory();
	
	/* load physics/geometric and LZerD/ZDOCK config, although it is not practically used anymore */
	this->decoy_program = options->get_decoy_program();
	this->score_type = options->get_score_type();

	// the "true" parameter supplied indicates that we should load the hydrogen-added version
	// of the pdbs. This is due to the fact that we need them for some of the scoring terms
	this->pdbs = load_pdbs(complex_name, chain_ids, directory, true);
	this->predictions = load_predictions(chain_ids, directory, decoy_program);

	/* number of generations and each generation's size */
	this->total_generations = options->get_generations() + this->current_generation;
	this->population_size = options->get_population();

	this->clash_threshold = options->get_clashes();
	this->stagnation_threshold = options->get_stagnation();
	this->stagnation_counter = stagnation_threshold;

	this->weights = options->get_weights();

	this->cluster_threshold = options->get_cluster_threshold();

	this->renewed_per_generation = options->get_renewed_amount();

	this->use_crossover = options->use_crossover();
	
	this->initial = generate_base_graph(chain_ids.size(),
                                      options->get_interactions());

  for(size_t edge_index = 0; edge_index < initial.E.size(); edge_index++) {
    md_edge edge = initial.E[edge_index];
    cout << edge.receptor << " " << edge.ligand << endl;  
  }

	// TODO: Similar to decoy_program and score_type this is not needed anymore
	// Find a way to remove this eventually (at least move to the main program)
	md_edge_set::set_scoring_scheme(score_type);

	/* if we are not resuming an execution, the population has not been initialized
	 * do it now that we've gathered all the configuration information */
	if(!options->resume_required())
	{
		population = generate_initial_population(predictions, initial.E, population_size);
	}
}

/*
 * This constructor is used when clustering an existing output file
 * We only initialize what's necessary for this purpose
 */
mdstate::mdstate(string input_filename, string directory, double cluster_cutoff)
{
	/* this method will initialize most of the instance variables needed from the file info */
	read(input_filename);
	// the "true" parameter supplied indicates that we should load the hydrogen-added version
	// of the pdbs. This is due to the fact that we need them for some of the scoring terms
	this->pdbs = load_pdbs(this->complex_name, this->chain_ids, directory, true);
	this->predictions = load_predictions(this->chain_ids, directory, this->decoy_program);

	this->cluster_threshold = cluster_cutoff;

	// TODO: Find a way to remove this eventually (at least move to the main program)
	md_edge_set::set_scoring_scheme(score_type);
}

/*
 * Augment the current population by applying crossover.
 * The parameter passed is used to output to the execution log
 */
void mdstate::crossover(ofstream& log_stream)
{
	log_stream << "Recombining..." << endl;
	vector<md_edge_set> recombined;
	// recombine pairs, in other words, 1 with 2, 3 with 4 and so on...
	for(size_t i = 0, j = 1; i < this->population.size() && j < this->population.size(); i+=2,j+=2)
	{
		if(!(population[i].equals(population[j]))) /* only recombine if both parents are different */
		{
			md_edge_set recombined_individual = recombine(population[i], population[j], this->predictions);
			recombined.push_back(recombined_individual);
		}
	}
	// join the recombined elements and the original population
	this->population.insert(this->population.end(), recombined.begin(), recombined.end());
	
	log_stream << "Finished Recombining..." << endl;
}

/*
 * Augment the current population by applying mutation.
 * The parameter passed is used to output to the execution log
 */
void mdstate::mutation(ofstream& log_stream)
{
	log_stream << "Mutating..." << endl;
	vector<md_edge_set> mutated;
	for(size_t k = 0; k < this->population.size(); k++)
	{
		md_edge_set individual = population[k];
		md_edge_set mutated_individual = mutate(individual, this->initial.E, this->predictions);
		mutated.push_back(mutated_individual);
	}
	// join the mutated elements and the original population
	this->population.insert(this->population.end(), mutated.begin(), mutated.end());

	log_stream << "Finished Mutating..." << endl;
}

/*
 * Score the current population and rank them accordingly. Note
 * that internally each graph's scoring method will not score
 * candidates that violate the clash threshold, however, they will
 * not be deleted yet.
 * The parameter passed is used to output to the execution log
 */
void mdstate::rank(ofstream& log_stream)
{
	log_stream << "Generating " << this->population.size() << " scores..." << endl;
	md_edge_set::update_score(&(this->population), this->pdbs, this->predictions, this->clash_threshold, this->weights);
	/*for(size_t w=0;w<this->population.size();w++)
	{
		population[w].update_score(this->pdbs, this->predictions, this->clash_threshold, this->weights);
	}*/
	log_stream << "Ranking..." << endl;
	sort(this->population.begin(), this->population.end(), md_edge_set::compare);
}

/*
 * Remove all the elements that violate the restrictions imposed
 * (i.e. clashes, duplicates and similarity identified by clustering)
 * The parameter passed is used to output to the execution log
 */
void mdstate::filter(ofstream& log_stream)
{
	log_stream << "Clash and duplicate prunning..." << endl;
	// after sorting prune those that have more clashes than what we want
	for(int w = this->population.size() - 1;
      (w >=0 && w < static_cast<int>(this->population.size())); w--)
	{
		if(population[w].get_clashes() > this->clash_threshold)
		{
			this->population.erase(this->population.begin() + w);
		}
		else if(w != (this->population.size() - 1))
		{ 
			/* i.e. if it's not the last,
			 * compare it to the next  and if they are the same don't add this one
			 * because we don't want duplicates */
			if(population[w].equals(population[w+1]))
			{
				this->population.erase(this->population.begin() + w);
			}
		}
	}

	// at this point the population only contains non-duplicate individuals 
	// and they meet the clash threshold.
	// Additionally, we'll cluster them if the cluster threshold is not zero
	if(this->cluster_threshold)
	{
		log_stream << "Clustering decoys..." << endl;
		this->cluster();
	}
}

/*
 * Takes the current clustering threshold and reduces the size of the population
 * according to the clustering algorithm implemented, keeping just the cluster centers
 *
 * return: the information corresponding to all cluster centers in case it needs to be
 * processed furthered (e.g. if we want to print the sizes of the clusters
 */
vector<cluster_neighborhood> mdstate::cluster()
{
    // Obtain the clusters. This static method will determine if MPI features can be used
    vector<cluster_neighborhood> centers = cluster_manager::cluster_population(&(this->population),
																&(this->pdbs), &(this->predictions),
																this->cluster_threshold);

	/* create a new population with the clustered elements that will replace the current population */
	vector<md_edge_set> clustered_population;
	for(size_t center_index = 0; center_index < centers.size(); center_index++)
	{
		clustered_population.push_back(population[centers[center_index].get_prediction_number()]);
	}

	this->population = clustered_population;

	/* the cluster centers are included in the new population starting from the one with the 
	 * highest number of neighbors, however, this does not necessarily represent the best scoring
	 * prediction. we sort in order to have the population correctl ranked
	 */
	sort(this->population.begin(), this->population.end(), md_edge_set::compare);
	
	return centers;
}

/*
 * Before finishing the current generation, make sure that we have enough
 * individuals to meet the minimum population size requirements. If too
 * many elements were filtered, then random new individuals are generated
 * although no processing will be done with them until the next iteration
 * The parameter passed is used to output to the execution log
 */
void mdstate::refill(ofstream& log_stream)
{
	if(this->population.size() < ((size_t)this->population_size))
	{
		// if we have at least one good individual reset the stagnation counter
		if(this->population.size() > 0)
		{
			this->stagnation_counter = this->stagnation_threshold;
			/*
			 * we want at least half of the population to come from the good candidates
			 * found. Thus we "clone" the good ones until we reach half of the population
			 * size expected. The rest will be filled with new individuals
			 */
			size_t good_individuals = this->population.size(); // the current size is the number of good candidates
			size_t current_duplicate_index = 0;
			while(this->population.size() < ((size_t)(this->population_size / 2)))
			{
				if(current_duplicate_index >= good_individuals)
				{
					// reset it
					current_duplicate_index = 0;
				}
				md_edge_set duplicate = population[current_duplicate_index];
				this->population.push_back(duplicate);
				current_duplicate_index++;
			}
		}
		else // we didn't have a good individual (at all). Decrease the stagnation counter
		{
			if(this->stagnation_counter-- == 0)
			{
				// Abort if we are over the limit
				log_stream << "Aborting. GA stagnation for " << this->stagnation_threshold
					<< " consecutive generations" << endl;
				log_stream.close();
				exit(1);
			}
		}
		/* we fill the rest up with totally new elements */
		size_t additional_individuals = this->population_size - this->population.size();
		log_stream << "Adding " << additional_individuals << " to complete the minimum population size..." << endl;
		this->renew(additional_individuals, true);
		/* since there are new elements we need to sort again */
		sort(this->population.begin(), this->population.end(), md_edge_set::compare);
	}
}

/*
 * This method adds N totally random new elements to the population
 */
void mdstate::renew(size_t additional_individuals, bool generate_score)
{
	vector<md_edge_set> additional_pop = generate_initial_population(this->predictions, this->initial.E, additional_individuals);
	/* update the score and then insert it */
	if(generate_score) { // score is computed only when it's a refill
		md_edge_set::update_score(&additional_pop, this->pdbs, this->predictions, this->clash_threshold, this->weights);
	}
	for(size_t additional_index = 0; additional_index < additional_pop.size(); additional_index++)
	{
		this->population.push_back(additional_pop[additional_index]);
	}
}

/*
 * Load the information contained in a MultiLZerD output file
 * and create the state that represents that file's contents
 */
void mdstate::read(string filename)
{
	// Read the information on the file
	ifstream ga_file_stream(filename.c_str(), ifstream::in);

	string decoy_program_str, score_type_str, chains_param, discard;
	// do it twice because we have "PDBID 1VCB", for example. The first one read PDBID to discard it
	ga_file_stream >> discard;
	ga_file_stream >> this->complex_name;
	// twice for the same reason
	ga_file_stream >> discard;
	ga_file_stream >> chains_param;

	// convert the raw comma-separated string to the vector representation
	this->chain_ids.clear();

	while(!chains_param.empty())
	{
		// get the following chain first
		size_t separator_pos = chains_param.find_first_of(',');
		string new_chain = chains_param.substr(0, separator_pos);
		this->chain_ids.push_back(new_chain);
		// and then remove the chain and the separator before the next iteration starts
		chains_param = separator_pos ==  string::npos ? "" : chains_param.substr(separator_pos + 1);
	}

	// Number of generations
	ga_file_stream >> discard;
	ga_file_stream >> this->current_generation;
	// Population size
	ga_file_stream >> discard;
	ga_file_stream >> this->population_size;
	// Clash threshold
	ga_file_stream >> discard;
	ga_file_stream >> this->clash_threshold;
	// Score type
	ga_file_stream >> discard;
	ga_file_stream >> score_type_str;
	this->score_type = score_type_str.compare(SCORETYPE_PAIRWISE) == 0 ?
					score_type_pairwise : score_type_physics;
	// Decoy program
	ga_file_stream >> discard;
	ga_file_stream >> decoy_program_str;
	this->decoy_program = decoy_program_str.compare(DECOYPROGRAM_LZERD) == 0 ?
					decoy_program_lzerd : decoy_program_zdock;

	this->population = load_population(ga_file_stream);
	// close the ga file because we don't need it anymore
	ga_file_stream.close();
}

/*
 * This method coordinates the different steps to execute one generation's
 * worth of simulation. In other words, it takes the current population and applies
 * the GA operations, ranking and filtering, and reflecting those changes by
 * modifying this instance's population
 *
 * log_stream: It receives a file stream that is used as a log to record the progress
 * throughout the different steps
 *
 * return: true if we still have generations left or false if we are finished (in this
 * case no actions are executed)
 */
bool mdstate::advance_generation(ofstream& log_stream)
{
	if(this->current_generation > this->total_generations)
	{
		return false;
	}

	time_t start_time, end_time;
	start_time = time(NULL);

	// else we still have generations left
	log_stream << "Generation " << this->current_generation  << " of " << this->total_generations << endl;
	if(this->current_generation > 0) {
		log_stream << "Previous generation executed in " << (this->previous_time / 60.0) << " minutes" << endl;
	}
	log_stream << "Stagnation counter " << this->stagnation_counter << endl;

	if(this->use_crossover)
	{
		this->crossover(log_stream);
	}

	this->mutation(log_stream);
	// add some totally random new elements
	if(this->renewed_per_generation) {
		log_stream << "Adding " << renewed_per_generation << " random individuals" << endl;
		this->renew(renewed_per_generation, false);
	}
	this->rank(log_stream);
	this->filter(log_stream);
	this->refill(log_stream);

	// only leave the number of individuals that we want to keep after each iteration
	this->population.erase(this->population.begin() + this->population_size, this->population.end());

	// finally update the generation counter and return true indicating that
	// the generation was processed correctly
	this->current_generation++;

	// update the time it took to execute this generation
	end_time = time(NULL);
	this->previous_time = difftime(end_time, start_time);

	return true;
}

/*
 * Write a MultiLZerD output file that describes the state of the simulation up to
 * this point. It can be used to generate the final status of the simulation or
 * to have detailed information about what happens after each generation
 */
void mdstate::write(ofstream& output_stream)
{
	string score_str = score_type == score_type_pairwise ? SCORETYPE_PAIRWISE : SCORETYPE_PHYSICS;
	string decoy_str = decoy_program == decoy_program_lzerd ? DECOYPROGRAM_LZERD : DECOYPROGRAM_ZDOCK;

	// Print a header describing the output
	output_stream << "PDBID\t" << this->complex_name << endl
			<< "Chains\t";
	// Print each of the chains followed by a comma
	for(size_t chain_index = 0; chain_index < this->chain_ids.size(); chain_index++)
	{
		output_stream << chain_ids[chain_index];
		if(chain_index < (this->chain_ids.size() - 1))
		{
			/* i.e. print a comma except when it's the last one */
			output_stream << ",";
		}
	}

	output_stream << endl << "Generations\t" << this->current_generation << endl
			<< "PopulationSize\t" << this->population_size << endl
			<< "ClashThreshold\t" << this->clash_threshold << endl
			<< "ScoreType\t" << score_str  << endl
			<< "DecoyProgram\t" << decoy_str << endl << endl;
	
	// Print the information corresponding to each individual
	for(size_t out = 0; out < this->population.size(); out++)
	{
		population[out].print(output_stream);
	}
}
