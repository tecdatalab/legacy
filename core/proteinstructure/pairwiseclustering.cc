#include <iostream>
#include <fstream>
#include <vector>
#include "clustering_options.h"
#include "cluster_neighborhood.h"
#include "cluster_manager.h"
#include "rmsd.h"
#include "lzerd_transformations.h"
#include "pdb.h"

using namespace std;

int main(int argc, char** argv)
{
	clustering_options options(argc,argv);
	if(!options.parse_successful())
	{
		cerr << "Usage: lzerd_cluster -R <receptor pdb file> -L <ligand pdb file> --input <LZerD output file> --cutoff <cluster distance>"
		<< " [--lrmsd] [--printsizes <cluster sizes output file>]" << endl;
		exit(EXIT_FAILURE);
	}

	string input_file = options.get_input_file();
	string ligand_file = options.get_ligand_file();
	string receptor_file = options.get_receptor_file();
	string sizes_file = options.get_sizes_file();
	bool print_sizes = !sizes_file.empty();
	double cutoff = options.get_cutoff();
	bool use_lrmsd = options.get_lrmsd_flag();

	// Load the information
	pdb receptor, ligand;

	// true, because we only want to load c-alpha atoms (all the others are ignored)
	read_protein(ligand_file, ligand, true);
	// false in this case because we use it to determine the LRMSD, which takes into account distances to any atom
	read_protein(receptor_file, receptor, false);

	// If the clustering should be done using LRMSD then we filter out all atoms that
	// are more than 10 angstroms away (from a receptor atom)
	if(use_lrmsd)
	{
		ligand.remove_distant_atoms(&(receptor.atoms), LIGAND_RMSD_THRESHOLD);
	}

	lzerd_transformations predictions(input_file);
	size_t total = predictions.get_size();

	// This collection will hold the transformed atoms (one position for each transformation)
	// to avoid recomputing the transformations
	vector<vector<atom> > transformed_ligands;

	vector<atom> transformed_ligand;
	vector<atom> base_ligand;

	for(size_t prediction_index = 0; prediction_index < total; prediction_index++)
	{
		transformed_ligand.clear();
		predictions.transform_atoms(ligand.atoms, transformed_ligand, prediction_index);
		transformed_ligands.push_back(transformed_ligand);
	}

	// Do the actual clustering
	cluster_manager cluster_handler(&transformed_ligands);
	cluster_handler.cluster(cutoff);
	// true means that we will use the highest ranked element of a cluster as its "center"
	vector<cluster_neighborhood> distances = cluster_handler.get_cluster_centers(true);

	// In case we need to output the sizes, we create an instance of a file output stream
	ofstream sizesout;
	if(print_sizes)
	{
		// delete it first, then open it
		remove(sizes_file.c_str());
		sizesout.open(sizes_file.c_str());
	}

	for(size_t distance_index = 0; distance_index < distances.size(); distance_index++)
	{
		predictions.print(cout, distances[distance_index].get_prediction_number());
		if(print_sizes)
		{
			// Print the size of this cluster as a single line (use the file specified in the parameters)
			sizesout << distances[distance_index].get_neighbor_count() << endl;
		}
	}

	// close the sizes file in case we used it
	if(print_sizes)
	{
		sizesout.close();
	}

	exit(EXIT_SUCCESS);
}
