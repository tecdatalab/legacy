#ifndef _LZERD_TRANSFORMATIONS_H_
#define _LZERD_TRANSFORMATIONS_H_

using namespace std;

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "transformations.h"

/*

#include "zdock.h"
#include "zdock_transformations.h"
#include <cstdlib>
*/

// Matrix that integrates the transformation info used by LZerD
typedef struct {
	double values[4][3];
} T43Matrix;

/*
 * This class inherits from the general interface "transformations"
 * and overrides the transformation method using LZerD-based decoys
 */
class lzerd_transformations : public transformations
{
	private:
		// Collection of predictions extracted from a LZerD output file. Each prediction is
		// stored as a matrix that integrates the transformation info
		vector<T43Matrix> predictions;
		// Each index represents the score of the n-th prediction in the LZerD output file
		vector<double> scores;
		// Base ligand rotation and translation. This is a random displacement that is optional
		// Note that we will not actually use this for multiple docking to avoid inconsistencies in
		// the pairwise predictions
		double rotL[3], transL[3];
		// These will be set to 1 if the LZerD file contains a random rotation and translation at the beginning
		char has_initial_transformation;

		/*
		 * Original version of a function that reads a LZerD output file
		 */
		void read_lzerd_data(string infile, vector<vector<double> >& Z, double* transL, double* rotL)
		{
			ifstream ifs;
			ifs.open(infile.c_str());
			if(!ifs)
			{
				cerr << "Could not open file: " << infile << endl;
				exit(EXIT_FAILURE);
			}

			string line;

			double t[13];
			while(getline(ifs, line))
			{
				if(line.find("LIG:") != string::npos)
				{
					istringstream iss(line);
					string s;
					iss >> s >> rotL[0] >> rotL[1] >> rotL[2] >> transL[0] >> transL[1] >> transL[2];
					has_initial_transformation = 1;
					continue;
				}

				vector<double> V;
				istringstream iss(line);
				for(int i=0;i<13;i++)
					iss >> t[i];
				for(int i=0;i<13;i++)
					V.push_back(t[i]);
				Z.push_back(V);
			}
			ifs.close();

			if(Z.size() == 0)
			{
				cerr << "Warning: no data to process in " << infile << endl;
			}
		}
		
		/*
		 * Original version of the method that applies the base rotation
		 */
		void apply_rotation(vector<atom>& atoms, double* rot)
		{
			double tx1, ty1, tz1, oldx, oldy, oldz;
			double rand1 = rot[0], rand2 = rot[1], rand3 = rot[2];
			vector<atom> NB;
			for(size_t i=0;i<atoms.size();i++)
			{
				atom RA = atoms[i];
				oldx = RA.axyz[0];
				oldy = RA.axyz[1];
				oldz = RA.axyz[2];

				rotate_coord(oldx, oldy, oldz, tx1, ty1, tz1, rand1, rand2, rand3);
				RA.axyz[0] = tx1;RA.axyz[1] = ty1;RA.axyz[2] = tz1;

				NB.push_back(RA);
			}
			atoms.clear();
			atoms = NB;
		}

		/*
		 * Original version of the method that applies the base translation
		 */
		void translate_atoms(vector<atom>& atoms, double* cog)
		{
			for(size_t i=0;i<atoms.size();i++)
			{
				for(int j=0;j<3;j++)
				{
					atoms[i].axyz[j] += cog[j];
				}
			}
		}

		/*
		 * Original version of the method used to rotate coordinates
		 */
		void rotate_coord(double oldX, double oldY, double oldZ,
					double& newX, double& newY, double& newZ,
					double phi, double theta, double psi)
		{
      double r11, r21, r31, r12, r22, r32, r13, r23, r33;

      r11 = cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi);
      r12 = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi);
      r13 = sin(psi)*sin(theta);
      r21 = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi);
      r22 = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi);
      r23 = cos(psi)*sin(theta);
      r31 = sin(theta)*sin(phi);
      r32 = -sin(theta)*cos(phi);
      r33 = cos(theta);

      newX = r11 * oldX + r12 * oldY + r13 * oldZ;
      newY = r21 * oldX + r22 * oldY + r23 * oldZ;
      newZ = r31 * oldX + r32 * oldY + r33 * oldZ;
		}

		/*
		 * Original version of the method used to apply a matrix transformation to
		 * a complete set of atoms
		 */
		void transform(T43Matrix R, vector<atom>& ain, vector<atom>& atrans)
		{
			for(size_t i=0;i<ain.size();i++)
			{
				atom t = ain[i];
				point_transform(t.axyz, R);
				atrans.push_back(t);
			}
		}

		/*
		 * Original version of the method used to transform a single point
		 * based on a matrix
		 */
		void point_transform (double* p, T43Matrix m)
		{
			double x, y, z;
			if(inverted)
			{
				/* first translate */
				p[0] += m.values[3][0];
				p[1] += m.values[3][1];
				p[2] += m.values[3][2];
				/* then apply the transformation */
				x = m.values[0][0]*p[0] + m.values[0][1]*p[1] + m.values[0][2]*p[2];
				y = m.values[1][0]*p[0] + m.values[1][1]*p[1] + m.values[1][2]*p[2];
				z = m.values[2][0]*p[0] + m.values[2][1]*p[1] + m.values[2][2]*p[2];
			}
			else /* do the regular transformation */
			{
				x = m.values[0][0]*p[0] + m.values[0][1]*p[1] + m.values[0][2]*p[2] + m.values[3][0];
				y = m.values[1][0]*p[0] + m.values[1][1]*p[1] + m.values[1][2]*p[2] + m.values[3][1];
				z = m.values[2][0]*p[0] + m.values[2][1]*p[1] + m.values[2][2]*p[2] + m.values[3][2];
			}
			p[0] = x; p[1] = y; p[2] = z;
		}

	public:
		/*
		 * Creates a new collection of LZerD transformations based on the file provided
		 */
		lzerd_transformations(string prediction_filename)
		{
			// These are the raw values read from the file. Each row represents a transformation
			// along with its score. The first 12 values will be used to create the transformation
			// matrix and the last one represents the score
			vector<vector<double> > raw_predictions;

			// initialize the overall transformation in case they're not provided by the file
			rotL[0] = rotL[1] = rotL[2] = transL[0] = transL[1] = transL[2] = 0;
			has_initial_transformation = 0; // this should be set to 1, in case there actually is one

			// read the data. Store the base transformations in instance variables (transL and rotL)
			// and load the raw predictions in order to generate the matrices
			read_lzerd_data(prediction_filename, raw_predictions, transL, rotL);

			// parse the raw information to extract the matrices and scores
			for(size_t i=0;i<raw_predictions.size();i++)
			{
				T43Matrix M;
				size_t k = 0;
				for(size_t l=0;l<3;l++)
				{
					for(size_t j=0;j<3;j++)
					{
						M.values[l][j] = raw_predictions[i][k];
						k++;
					}
				}
				M.values[3][0] = raw_predictions[i][9];
				M.values[3][1] = raw_predictions[i][10];
				M.values[3][2] = raw_predictions[i][11];

				// store the matrix
				predictions.push_back(M);
				// and now extract the score
				scores.push_back(raw_predictions[i][12]);
			}
		}

		/*
		 * This method is overridden in order to implement atom transformations based on zdock decoys
		 * original represents the atoms before transformation
		 * transformed is assumed to be an empty vector where the transformed atoms will be added
		 * transformation_number represents the index of the transformation used,
		 * since this instances contains a collection of transformations
		 */
		void transform_atoms(vector<atom>& original, vector<atom>& transformed, size_t transformation_number)
		{
			// first make a copy because we need to apply the base translation and rotation
			vector<atom> copy = original;
			// apply base translation and rotation
			if(has_initial_transformation)
			{
				translate_atoms(copy, transL);
				apply_rotation(copy, rotL);
			}
			// and now get the transformation matrix and apply it
			transform(predictions[transformation_number], copy, transformed);
		}

		/*
		 * Returns the ZDOCK score assigned to the prediction number specified by the index
		 */
		double get_score(size_t transformation_number)
		{
			return scores[transformation_number];
		}

		/*
		 * Returns the number of predictions held by this instance
		 */
		size_t get_size()
		{
			return predictions.size();
		}

		/*
		 * Prints the n-th ranked prediction only (prediction_number)
		 * It outputs the transformation followed by the score, in the same
		 * format as the original LZerD output
		 */
		void print(ostream& output_stream, size_t prediction_number)
		{
			T43Matrix M = predictions[prediction_number];
			output_stream << fixed << left << setprecision(3);
			for(size_t l=0;l<3;l++)
			{
				for(size_t j=0;j<3;j++)
				{
					output_stream << setw(6) << M.values[l][j] << "  ";
				}
			}
			output_stream << setw(6) << M.values[3][0] << "  "
					<< setw(6) << M.values[3][1] << "  "
					<< setw(6) << M.values[3][2] << "  "
					<< setw(6) << scores[prediction_number]
					<< endl;
		}

		/*
		 * Prints all the information corresponding to each transformation
		 * to the output stream provided
		 */
		void print(ostream& output_stream)
		{
			for(size_t t = 0; t < predictions.size(); t++)
			{
				output_stream << "#" << t << ":\n";
				for(size_t i = 0; i < 4; i++)
				{
					for(size_t j = 0; j < 3; j++)
					{
						output_stream << predictions[t].values[i][j] << " ";
					}
					output_stream << endl;
				}
				output_stream << endl;
			}
		}

		/*
		 * Inverting a LZerD transformation means inverting the rotation matrix
		 * and changing the sign of the translation values. The initial random translation
		 * and rotation of the ligand are ignored because we don't use them when generating decoys
		 */
		void invert()
		{
			/* Invert all of the predictions */
			for(size_t prediction = 0; prediction < predictions.size(); prediction++)
			{
				/* Create a new matrix that will replace the current one */
				T43Matrix inverted;
				/* Get the current matrix to use it as a base */
				T43Matrix old = predictions[prediction];
				/* Copy the translation values and change their sign */
				inverted.values[3][0] = -old.values[3][0];
				inverted.values[3][1] = -old.values[3][1];
				inverted.values[3][2] = -old.values[3][2];

				/* Calculate the inverse of the 3x3 rotation portion of the matrix
				 * by using the determinant definition for 3x3 and 2x2 matrices */
				double det = (old.values[0][0]*(old.values[1][1]*old.values[2][2] - old.values[2][1]*old.values[1][2]) -
						old.values[0][1]*(old.values[1][0]*old.values[2][2] - old.values[2][0]*old.values[1][2]) +
						old.values[0][2]*(old.values[1][0]*old.values[2][1] - old.values[2][0]*old.values[1][1]));

				if(det == 0.0)
				{
					/* In general, this shouldn't happen, since a rotation matrix should
					 * have non-zero determinant, however, for chains that are not interacting this might be the case
					 * We just use the transpose instead of the standard definition*/
					inverted.values[0][0] = old.values[0][0];
					inverted.values[1][1] = old.values[1][1];
					inverted.values[2][2] = old.values[2][2];
					
					inverted.values[0][1] = old.values[1][0];
					inverted.values[0][2] = old.values[2][0];

					inverted.values[1][0] = old.values[0][1];
					inverted.values[1][2] = old.values[2][1];

					inverted.values[2][0] = old.values[0][2];
					inverted.values[2][1] = old.values[1][2];
				}
				else
				{
					inverted.values[0][0] = (old.values[1][1]*old.values[2][2] - old.values[1][2]*old.values[2][1])/det;
					inverted.values[0][1] = (old.values[2][1]*old.values[0][2] - old.values[0][1]*old.values[2][2])/det;
					inverted.values[0][2] = (old.values[0][1]*old.values[1][2] - old.values[1][1]*old.values[0][2])/det;
					inverted.values[1][0] = (old.values[1][2]*old.values[2][0] - old.values[1][0]*old.values[2][2])/det;
					inverted.values[1][1] = (old.values[0][0]*old.values[2][2] - old.values[2][0]*old.values[0][2])/det;
					inverted.values[1][2] = (old.values[1][0]*old.values[0][2] - old.values[0][0]*old.values[1][2])/det;
					inverted.values[2][0] = (old.values[1][0]*old.values[2][1] - old.values[2][0]*old.values[1][1])/det;
					inverted.values[2][1] = (old.values[2][0]*old.values[0][1] - old.values[0][0]*old.values[2][1])/det;
					inverted.values[2][2] = (old.values[0][0]*old.values[1][1] - old.values[0][1]*old.values[1][0])/det;
				}

				/* replace the old matrix by the inverted one */
				predictions[prediction] = inverted;
			}
			/* change the inverted flag */
			inverted = !inverted;
		}

	};
#endif
