#ifndef _MDSTATE_H_
#define _MDSTATE_H_

#include <fstream>
#include <vector>
#include <ctime>
#include "mdga_options.h"
#include "mdgraph.h"
#include "pdb.h"
#include "transformations.h"
#include "cluster_neighborhood.h"
#include "cluster_manager.h"

/*
 * An instance of this class represents the state of a multiple docking simulation
 * after a certain number of generations have been executed.
 * This same class provides the necessary methods to initialize a new simulation
 * and to advance from the current state to the next generation, by applying the
 * GA operations defined
 */
class mdstate
{
	private:
		/*
		 * This instance variable stores the graphs before an iteration is
		 * started and contains the final elements kept for the next iteration
		 * after applying all the different genetic algorithm operations
		 */
		vector<md_edge_set> population;
		/* 
		 * This number is different from population.size(); this int value represents
		 * the number of individuals that must be kept after each generation
		 */
		int population_size;
		/*
		 * This is the number of totally new elements generated per generation
		 * The intention behind this is to improve the sampling
		 */
		int renewed_per_generation;
		/* These two keep track of the generation index and the last generation */
		int current_generation;
		int total_generations;
		/* 
		 * The complex name and chain ID's are used to write the header of the output file
		 * that reflects the status of this mdstate
		 */
		string complex_name;
		vector<string> chain_ids;
		/* Holds the input atomic structures that will be re-assembled */
		vector<pdb> pdbs;
		/* 
		 * Stores the translation/rotation transformation corresponding to all pairwise predictions generated
		 * by LZerD
		 */
		vector< vector<transformations*> > predictions;
		/*
		 * Elements with more than this amount of clashes will be deleted from population at the end of
		 * each iteration
		 */
		int clash_threshold;
		/*
		 * If no suitable candidates are found after this number of iterations, the simulation stops
		 */
		int stagnation_threshold;
		/* 
		 * stagnation_threshold is a read-only value, this counter is the variable that is actually
		 * used during the execution
		 */
		int stagnation_counter;
		/*
		 * A linear combination of 12 terms is used to assess the fitness of each of the individuals in the
		 * population. Each execution of the multiple docking program can use a different set of weights.
		 * This pointer corresponds to the 12 weights used to generate the ranking score
		 */
		double* weights;
		/* 
		 * These 2 configuration values were used at the beginning of the development, however,
		 * only decoy_program=LZerD and score_type=physics are used at this point.
		 * They will be left in the code in case there is a change in the future
		 */
		decoy_program_t decoy_program;
		score_type_t score_type;
		/*
		 * Clustering threshold that takes a zero value if no clustering should be performed
		 * or the distance threshold, in terms of C-alpha RMSD, to group decoys
		 */
		double cluster_threshold;
		/* The use of crossover is optional (if true it will be used) */
		bool use_crossover;
		/* 
		 * Store the time it took the last generation to finish. This time will be output to the
		 * log whenever a new generation starts
		 */
		double previous_time;
		/* 
		 * This graph is used to generate new random individuals. In the standard configuration
		 * this should be a fully connected graph, however, if restrictions need to be imposed it would
		 * be as easy as deleting some of the edges from the initial graph, and none of the individuals
		 * will have them
		 */
		md_graph initial;

		/*
		 * Augment the current population by applying crossover.
		 * The parameter passed is used to output to the execution log
		 */
		void crossover(ofstream& log_stream);
		/*
		 * Augment the current population by applying mutation.
		 * The parameter passed is used to output to the execution log
		 */
		void mutation(ofstream& log_stream);
		/*
		 * Score the current population and rank them accordingly. Note
		 * that internally each graph's scoring method will not score
		 * candidates that violate the clash threshold, however, they will
		 * not be deleted yet.
		 * The parameter passed is used to output to the execution log
		 */
		void rank(ofstream& log_stream);
		/*
		 * Remove all the elements that violate the restrictions imposed
		 * (i.e. clashes, duplicates and similarity identified by clustering)
		 * The parameter passed is used to output to the execution log
		 */
		void filter(ofstream& log_stream);
		/*
		 * Before finishing the current generation, make sure that we have enough
		 * individuals to meet the minimum population size requirements. If too
		 * many elements were filtered, then random new individuals are generated
		 * although no processing will be done with them until the next iteration
		 * The parameter passed is used to output to the execution log
		 */
		void refill(ofstream& log_stream);

		/*
		 * This method adds N totally random new elements to the population
		 */
		void renew(size_t additional_individuals, bool generate_score);
		/*
		 * Load the information contained in a MultiLZerD output file
		 * and create the state that represents that file's contents
		 */
		void read(string filename);
	public:
		/*
		 * Initializes the simulation state based on the configuration provided
		 * Note that this can be either the start of a new simulation or we can
		 * be resuming a previous one
		 */
		mdstate(mdga_options* options);
		/*
		 * This constructor is used when clustering an existing output file
		 * We only initialize what's necessary for this purpose
		 */
		mdstate(string input_file, string directory, double cluster_cutoff);
		/*
		 * Returns the number of generations that have been executed so far
		 */
		int get_current_generation()
		{
			return this->current_generation;
		}
		/*
		 * Returns a reference to the PDB's loaded in the constructor
		 */
		vector<pdb>* get_pdbs()
		{
			return &pdbs;
		}
		/*
		 * Returns a reference to the transformations loaded in the constructor
		 */
		vector< vector<transformations*> >* get_predictions()
		{
			return &predictions;
		}
		/*
		 * Returns the value of the clash cutoff used in the simulation
		 */
		int get_clash_cutoff()
		{
			return clash_threshold;
		}
		/*
		 * Returns a reference to the 12-term vector used
		 * to generate the linear combination of scores
		 */
		double* get_weights()
		{
			return weights;
		}
		/*
		 * Number of elements maintained as current population
		 */
		int get_population_size()
		{
			return population_size;
		}
		/*
		 * This method coordinates the different steps to execute one generation's
		 * worth of simulation. In other words, it takes the current population and applies
		 * the GA operations, ranking and filtering, and reflecting those changes by
		 * modifying this instance's population
		 * log_stream: It receives a file stream that is used as a log to record the progress
		 * throughout the different steps
		 * 
		 * return: true if we still have generations left or false if we are finished (in this
		 * case no actions are executed)
		 */
		bool advance_generation(ofstream& log_stream);
		/*
		 * Write a MultiLZerD output file that describes the state of the simulation up to
		 * this point. It can be used to generate the final status of the simulation or
		 * to have detailed information about what happens after each generation
		 */
		void write(ofstream& output_stream);
		/*
		 * Takes the current clustering threshold and reduces the size of the population
		 * according to the clustering algorithm implemented, keeping just the cluster centers
		 * 
		 * return: the information corresponding to all cluster centers in case it needs to be
		 * processed furthered (e.g. if we want to print the sizes of the clusters
		 */
		vector<cluster_neighborhood> cluster();
};

#endif
