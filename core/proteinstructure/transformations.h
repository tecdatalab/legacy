#ifndef _TRANSFORMATIONS_H_
#define _TRANSFORMATIONS_H_
#include <vector>
#include <iostream>
#include "atom.h"

/*
 * This is an abstract class that will serve as base for different implementations
 * of transformations, using decoy generation programs like ZDOCK and LZerD.
 * Since each of those programs can potentially have different ways of expressing the transformations
 * the implementation can vary considerably.
 * Each instance of this class should contain a collection of rotation and translation parameters, each of which
 * represents a decoy.
 * Each subclass must provide a means of transforming atoms using any of the transformations stored in this
 * instance and also a way of retrieving the score of pairwise decoys. Additionally, there must be a method to
 * retrieve the number of predictions stored by a transformations instance.
 * For practical purposes, there should be one instance of this class for each ZDOCK or LZerD output file,
 * representing pairwise decoys of 2 structures.
 */
class transformations
{
	public:
		bool inverted;
	public:
		transformations()
		{
			/* initially this instance does not represent an inverted transformation */
			inverted = false;
		}

		/* Added only because moffet uses an old compiler that throws a warning if the destructor is not present */
		virtual ~transformations()
		{
		}

		/*
		 * This method needs to be overridden by child classes.
		 * original represents the atoms before transformation
		 * transformed is assumed to be an empty vector where the transformed atoms will be added
		 * transformation_number represents the index of the transformation used,
		 * since this instances contains a collection of transformations
		 * The inverted boolean value instructs the method to translate first and then rotate (normally it would
		 * be rotate then translate). This happens when we are using one of the inverted transformations
		 */
		virtual void transform_atoms(vector<atom>& original, vector<atom>& transformed, size_t transformation_number) = 0;
		/*
		 * This method should return the numeric score assigned to a pairwise combination of receptor and the ligand after applying
		 * the transformation indexed by this parameter
		 */
		virtual double get_score(size_t transformation_number) = 0;
		/*
		 * Return the total number of predictions
		 */
		virtual size_t get_size() = 0;
		/*
		 * Prints all the information corresponding to each transformation
		 * to the output stream provided
		 */
		virtual void print(ostream& output_stream) = 0;
		/*
		 * Inverts the translation and transformation information
		 * of this instance. The idea behind this is that normally decoys are generated for
		 * 2 chains in one direction, for example, setting chain A as receptor and B as ligand.
		 * However, it is possible that a graph has an edge with direction B->A. To apply that
		 * transformation we can take the one corresponding to A->B and invert it to get the desired
		 * result
		 */
		virtual void invert() = 0;
};
#endif
