#ifndef _MDGA_OPTIONS_H_
#define _MDGA_OPTIONS_H_

#include <cstring>
#include <getopt.h>

#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include "md.h"
#include "mdgraph.h"
#include "score_weights_options.h"

using namespace std;

#define GENERATIONS_DEFAULT "1000"
#define POPULATION_DEFAULT "200"
#define CLASHES_DEFAULT "250"
#define STAGNATION_DEFAULT "5"
#define DIRECTORY_DEFAULT ""

#define OUTPUT_FILE_SUFFIX ".ga.out"
#define LOG_FILE_SUFFIX ".log"

#define PDBID_OPTION "pdbid"
#define PDBID_CHAR 'p'
#define CHAINS_OPTION "chains"
#define CHAINS_CHAR 'c'
#define GENERATIONS_OPTION "generations"
#define GENERATIONS_CHAR 'g'
#define POPULATION_OPTION "population"
#define POPULATION_CHAR 's'
#define CLASHES_OPTION "clashes"
#define CLASHES_CHAR 'l'
#define STAGNATION_OPTION "stagnation"
#define STAGNATION_CHAR 'a'
#define OUTPUT_PREFIX_CHAR 'o'
#define DIRECTORY_CHAR 'd'

#define SCORETYPE_OPTION "scoretype"
#define SCORETYPE_CHAR 'r'

#define DECOYPROGRAM_OPTION "decoyprogram"
#define DECOYPROGRAM_CHAR 'y'

#define RESUME_OPTION "resume"
#define RESUME_CHAR 'u'

#define DETAILED_OPTION "detailed"
#define DETAILED_CHAR 'e'

#define CROSSOVER_OPTION "crossover"
#define CROSSOVER_CHAR 'v'

#define CLUSTER_OPTION "cluster"
#define CLUSTER_CHAR 't'

#define RENEW_OPTION "renew"
#define RENEW_CHAR 'n'

#define INTERACTIONS_OPTION "interactions"
#define INTERACTIONS_CHAR 'i'

// o is the prefix used for the output files and d is the input directory, where the pdb and prediction files are located
#define SHORT_OPTIONS "o:d:"

// Option 1:
//--pdbid <PDB ID> --chains <A,B,C,etc> --generations <repetitions> --population <size>
//	--clashes <num> -o <output_prefix> -d <input_directory>
// Option 2:
//-resume <oldoutput.ga.out> --generations <repetitions> --population <size>
//	--clashes <num> -o <output_prefix> -d <input_directory>
static struct option mdga_long_options[] = 
{
	{PDBID_OPTION, required_argument, 0, PDBID_CHAR},
	{CHAINS_OPTION, required_argument, 0, CHAINS_CHAR},
	{GENERATIONS_OPTION, required_argument, 0, GENERATIONS_CHAR},
	{POPULATION_OPTION, required_argument, 0, POPULATION_CHAR},
	{CLASHES_OPTION, required_argument, 0, CLASHES_CHAR},
	{SCORETYPE_OPTION, required_argument, 0, SCORETYPE_CHAR},
	{DECOYPROGRAM_OPTION, required_argument, 0, DECOYPROGRAM_CHAR},
	{STAGNATION_OPTION, required_argument, 0, STAGNATION_CHAR},
	{RESUME_OPTION, required_argument, 0, RESUME_CHAR},
	{DETAILED_OPTION, no_argument, 0, DETAILED_CHAR},
	{CROSSOVER_OPTION, no_argument, 0, CROSSOVER_CHAR},
	{WEIGHTS_OPTION, required_argument, 0, WEIGHTS_CHAR},
	{CLUSTER_OPTION, required_argument, 0, CLUSTER_CHAR},
	{RENEW_OPTION, required_argument, 0, RENEW_CHAR},
	{INTERACTIONS_OPTION, required_argument, 0, INTERACTIONS_CHAR},
	{0,0,0,0}
};


// This class parses the command line arguments supplied to the multiple docking ga program
class mdga_options : public score_weights_options
{
	public:
		// Only constructor in this class. It uses the input arguments in order to call getopt_long
		// to help the parsing process. After the instance is created, the methods provided by this class
		// will tell the caller if the process was successful and what is the value of the different configuration
		// parameters used
		mdga_options(int argc, char** argv) : score_weights_options()
		{
			// Initialize argument values to default (if they're optional) or to "" (if they're required)
			main_complex_name = chains_arg = output_file_prefix = resume_param = cluster_param = renew_param = interactions = "";
			clashes_param = CLASHES_DEFAULT;
			directory = DIRECTORY_DEFAULT;
			generations_param = GENERATIONS_DEFAULT;
			population_size_param = POPULATION_DEFAULT;
			score_type = SCORETYPE_PHYSICS;
			decoy_program = DECOYPROGRAM_LZERD;
			stagnation_param = STAGNATION_DEFAULT;
			// by default we don't provide detailed output, resume is disabled and we don't crossover
			detailed = false;
			resume = false;
			crossover = false;
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c; // c represents the character associated with an option
				c = getopt_long(argc, argv, SHORT_OPTIONS, mdga_long_options, &option_index);
				if(c == -1)
				{
					options_left = false;
				}
				else
				{
					// determine what type of option we got
					switch(c)
					{
						case PDBID_CHAR:
							main_complex_name = optarg;
							break;
						case CHAINS_CHAR:
							chains_arg = optarg;
							break;
						case GENERATIONS_CHAR:
							generations_param = optarg;
							break;
						case POPULATION_CHAR:
							population_size_param = optarg;
							break;
						case CLASHES_CHAR:
							clashes_param = optarg;
							break;
						case OUTPUT_PREFIX_CHAR:
							output_file_prefix = optarg;
							break;
						case DIRECTORY_CHAR:
							directory = optarg;
							break;
						case SCORETYPE_CHAR:
							score_type = optarg;
							// check that it's one of the allowed types
							if(score_type.compare(SCORETYPE_PAIRWISE) != 0 &&
								score_type.compare(SCORETYPE_PHYSICS) != 0)
							{
								// it's not equal to the allowed options. Halt.
								cerr << "Allowed values for --" << SCORETYPE_OPTION << " are "
									<< SCORETYPE_PAIRWISE << " and " << SCORETYPE_PHYSICS << endl;
								options_left = false;
								parse_error = true;
							}
							break;
						case DECOYPROGRAM_CHAR:
							decoy_program = optarg;
							// check that it's one of the allowed types
							if(decoy_program.compare(DECOYPROGRAM_LZERD) != 0 &&
								decoy_program.compare(DECOYPROGRAM_ZDOCK) != 0)
							{
								// it's not equal to the allowed options. Halt.
								cerr << "Allowed values for --" << DECOYPROGRAM_OPTION << " are "
									<< DECOYPROGRAM_LZERD << " and " << DECOYPROGRAM_ZDOCK << endl;
								options_left = false;
								parse_error = true;
							}
							break;
						case STAGNATION_CHAR:
							stagnation_param = optarg;
							break;
						case RESUME_CHAR:
							resume_param = optarg;
							resume = true;
							break;
						case DETAILED_CHAR:
							detailed = true;
							break;
						case CROSSOVER_CHAR:
							crossover = true;
							break;
						case WEIGHTS_CHAR:
							this->set_weights(optarg);
                                                        break;
						case CLUSTER_CHAR:
							cluster_param = optarg;
							break;
						case RENEW_CHAR:
							renew_param = optarg;
							break;
            case INTERACTIONS_CHAR:
              interactions = optarg;
              break;
						case '?': //option not recognized (an error is automatically printed)
						default:
							options_left = false;
							parse_error = true;
							break;
					}
				}

			}
		}

		// the parse was successful if the errors flag is false and if all required arguments were supplied
		bool parse_successful()
		{
			return !parse_error && !output_file_prefix.empty() && this->get_weights() &&
				((!main_complex_name.empty() && !chains_arg.empty()) ||
					resume
				);
		}

		// Returns the PDBID that identifies the complex being processed
		string get_main_complex_name()
		{
			return main_complex_name;
		}

		// Parses the --chains option and returns a vector with the different chain identifiers
		vector<string> get_chain_ids()
		{
			string chains_param = chains_arg;
			vector<string> chain_ids;
			while(!chains_param.empty())
			{
				// get the following chain first
				size_t separator_pos = chains_param.find_first_of(',');
				string new_chain = chains_param.substr(0, separator_pos);
				chain_ids.push_back(new_chain);
				// and then remove the chain and the separator before the next iteration starts
				chains_param = separator_pos ==  string::npos ? "" : chains_param.substr(separator_pos + 1);
			}
			return chain_ids;
		}

		// Get a single string representation of the chains
		string get_all_chains()
		{
			return chains_arg;
		}

		// returns the number of generations provided or its default value, if it wasn't provided
		int get_generations()
		{
			return atoi(generations_param.c_str());
		}

		// returns the population size provided or its default value, if it wasn't provided
		int get_population()
		{
			return atoi(population_size_param.c_str());
		}
		// returns the clash threshold or its default value, if it wasn't provided
		int get_clashes()
		{
			return atoi(clashes_param.c_str());
		}
		// returns the stagnation threshold or its default value, if it wasn't provided
		int get_stagnation()
		{
			return atoi(stagnation_param.c_str());
		}

		// returns the name of the input directory provided or its default value, if it wasn't provided
		string get_input_directory()
		{
			return directory;
		}

		// Returns the name that should be used to create the output file, that contains the decoys specifications
		string get_output_filename()
		{
			return output_file_prefix + OUTPUT_FILE_SUFFIX;
		}

		// Returns the name of the log file that should be used, based on the prefix provided in the parameters
		string get_log_filename()
		{
			return output_file_prefix + LOG_FILE_SUFFIX;
		}

		// Returns the type of scoring used
		score_type_t get_score_type()
		{
			if(score_type.compare(SCORETYPE_PAIRWISE) == 0)
			{
				return score_type_pairwise;
			}
			else // SCORETYPE_PHYSICS
			{
				return score_type_physics;
			}
		}

		// Returns the option that establishes the program used to generate pairwise decoys
		decoy_program_t get_decoy_program()
		{
			if(decoy_program.compare(DECOYPROGRAM_LZERD) == 0)
			{
				return decoy_program_lzerd;
			}
			else // DECOYPROGRAM_ZDOCK
			{
				return decoy_program_zdock;
			}
		}

		// If provided the GA should restart execution using the population in this file
		string get_resumed_filename()
		{
			return resume_param;
		}

		// True is a file was provided as initial population, possibly the output from a previous execution
		bool resume_required()
		{
			return resume;
		}

		// Returns true if the output needs to be detailed for each generation
		bool is_detailed()
		{
			return detailed;
		}
		// Returns true if we need to do crossover
		bool use_crossover()
		{
			return crossover;
		}
		// Returns zero if no threshold was set or the value specified for clustering
		double get_cluster_threshold()
		{
			return cluster_param.empty() ? 0 : atof(cluster_param.c_str());
		}
		// Returns zero if no threshold was set or the value specified for clustering
		int get_renewed_amount()
		{
			return renew_param.empty() ? 0 : atoi(renew_param.c_str());
		}

    // Deals with the --interactions param and returns a series of paired
    // indices that represent which pairwise interactions are allowed.
    // For example, if --chains "G,B,J,M" and --interactions "B-M,G-B,G-M,B-J"
    // the result contains (1,3) (0,1) (0,3) (1,2) to reflect
    // the indices of the chain letters.
    vector<pair<md_vertex,md_vertex> > get_interactions()
    {
      vector<pair<md_vertex,md_vertex> > interacting_pairs;
      if(!interactions.empty()) {
        // Using the information from --chains create a tmp hash_map
        // to quickly determine the indexx for each chain letter
        map<string,md_vertex> chains_map;
        vector<string> chains = this->get_chain_ids();
        for(size_t chain_index = 0; chain_index < chains.size(); chain_index++) {
          chains_map[chains[chain_index]] = static_cast<md_vertex>(chain_index);
        }
        // Copy to be modifed by the tokenizer
        char* tokenized = new char[interactions.size() + 1];
        strncpy(tokenized, interactions.c_str(), interactions.size());
        tokenized[interactions.size()] = 0;
        char* current_chain = strtok(tokenized, ",-");
        string first_chain, second_chain;
        while(current_chain != NULL) {
          first_chain = current_chain;
          // move to the second chain
          current_chain = strtok(NULL, ",-");
          second_chain = current_chain;
          md_vertex first_vertex = chains_map[first_chain];
          md_vertex second_vertex = chains_map[second_chain];

          interacting_pairs.push_back(make_pair(first_vertex, second_vertex));
          // look for the next first part of a pair
          current_chain = strtok(NULL, ",-");
        }
        delete [] tokenized;
      }
      return interacting_pairs;
    }

	private:
		// These variables hold the parameters that come directly from the command line further processing is done by the methods
		string main_complex_name;
		string chains_arg;
		string generations_param;
		string population_size_param;
		string clashes_param;
		string output_file_prefix;
		string log_filename;
		string output_filename;
		string directory;
		string score_type;
		string decoy_program;
		string stagnation_param;
		string resume_param;
		string cluster_param;
		string renew_param;
    string interactions;
		bool detailed;
		bool resume;
		bool crossover;
		// set to false in case that unrecognized options or errors are found by getopt_long
		bool parse_error;
};

#endif
