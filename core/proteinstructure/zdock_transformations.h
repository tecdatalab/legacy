#ifndef _ZDOCK_TRANSFORMATIONS_H_
#define _ZDOCK_TRANSFORMATIONS_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include "transformations.h"
#include "zdock.h"

/*
using namespace std;

#include "zdock.h"
#include "zdock_transformations.h"
#include <cmath>
#include <cstdlib>
*/


/*
 * This class inherits from the general interface "transformations"
 * and overrides the transformation method using zdock-based decoys
 */
class zdock_transformations : public transformations
{
	private:
		// Collection of predictions extracted from a zdock output file
		zdata predictions;

		/*
		 * Original version of a function that reads a zdock output file and creates a zdata instance
		 */
		void read_zdock_data(string infile, string& fixed, string& mov, zdata& Z)
		{
			ifstream ifs;
			ifs.open(infile.c_str());
			if(!ifs)
			{
				cerr << "Could not open file: " << infile << endl;
				exit(EXIT_FAILURE);
			}
    

			string s1, s2;
    
			string line;
			int lnc = 0;
			double t1[3];
			int t2[3];
			double score;
			while(getline(ifs, line))
			{
				lnc = lnc + 1;
				istringstream iss(line);
				if(lnc == 1)
				{
					iss >> Z.N >> Z.SPACING;
				}
				if(lnc == 2)
				{
					iss >> Z.RANDVECT[0] >> Z.RANDVECT[1] >> Z.RANDVECT[2];
				}
				if(lnc == 3)
				{
					iss >> s1 >> Z.RVECT[0] >> Z.RVECT[1] >> Z.RVECT[2];
					//fixed = get_file_name(s1) + ".pqr";
					fixed = s1;
				}
				if(lnc == 4)
				{
					iss >> s2 >> Z.LVECT[0] >> Z.LVECT[1] >> Z.LVECT[2];
					//mov = get_file_name(s2) + ".pqr";
					mov = s2;
				}
				if(lnc > 4)
				{
            				iss >> t1[0] >> t1[1] >> t1[2] >> t2[0] >> t2[1] >> t2[2] >> score;
					zdock nz(t1, t2, score);
					Z.ZD.push_back(nz);
				}
			}
			ifs.close();

			if(Z.ZD.size() == 0)
			{
				cerr << "Error: no data to process" << endl;
				exit(EXIT_FAILURE);
			}
		}

		/**********************************************/
		/* Function rotateAtom  */
		/*  rotates around 3 euler angles  */
		/**********************************************/
		void rotateAtom (double oldX, double oldY, double oldZ,
                		 double& newX, double& newY, double& newZ,
		                 double psi, double theta, double phi)
		{
			float r11, r21, r31, r12, r22, r32, r13, r23, r33;
			
			r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi);
			r21 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi);
			r31 = sin(theta)*sin(phi);

			r12 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi);
			r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi);
			r32 = sin(theta)*cos(phi);

			r13 = sin(psi)*sin(theta);
			r23 = -cos(psi)*sin(theta);
			r33 = cos(theta);

			newX = r11 * oldX + r12 * oldY + r13 * oldZ;
			newY = r21 * oldX + r22 * oldY + r23 * oldZ;
			newZ = r31 * oldX + r32 * oldY + r33 * oldZ;
		}


		/*
		 * This is the original version of the function used in pairwise docking. It was left as is
		 * in case a future comparison needs to be made (against the old code). Only transform_atoms should be invoked
		 * from outside the class
		 */
		void transform_coords(vector<atom>& A, vector<atom>& atrans, double rand1,
			double rand2, double rand3, double r1, double r2, double r3,
			double l1, double l2, double l3, double a1, double a2, double a3,
			int t1, int t2, int t3, int N, double spacing)
		{
			//cerr << "A=" << A.size() << endl;
			double tx1, ty1, tz1, tx2, ty2, tz2, oldx, oldy, oldz;
			for(size_t i=0;i<A.size();i++)
			{
				atom RA = A[i];
				oldx = RA.axyz[0] - l1;
				oldy = RA.axyz[1] - l2;
				oldz = RA.axyz[2] - l3;

				rotateAtom(oldx, oldy, oldz, tx1, ty1, tz1, rand1, rand2, rand3);
				rotateAtom(tx1, ty1, tz1, tx2, ty2, tz2, a1, a2, a3);

				/* adjust so coordinates are in the box */
				if (t1 >= N/2) t1 -= N;
				if (t2 >= N/2) t2 -= N;  
				if (t3 >= N/2) t3 -= N;

				//grid-adjusted coords
				RA.axyz[0] = tx2 - t1*spacing + r1;
				RA.axyz[1] = ty2 - t2*spacing + r2;
				RA.axyz[2] = tz2 - t3*spacing + r3;

				atrans.push_back(RA);
			}
		}

	public:
		/*
		 * Creates a new collection of zdock transformations based on the file provided
		 */
		zdock_transformations(string prediction_filename)
		{
			// these are two strings output by the read_zdock_data function but we don't use them (however they are required)
			string fixed, mov;
			read_zdock_data(prediction_filename, fixed, mov, predictions);
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
			// TODO: fix this method and transform_coords to use the "inverted" instance variable
			// At this point this has only been corrected for the LZerD case
			zdock prediction_transformation = predictions.ZD[transformation_number];

			transform_coords(original, transformed,
					predictions.RANDVECT[0], predictions.RANDVECT[1], predictions.RANDVECT[2],
					predictions.RVECT[0], predictions.RVECT[1], predictions.RVECT[2],
					predictions.LVECT[0], predictions.LVECT[1], predictions.LVECT[2],
					prediction_transformation.CANGLE[0], prediction_transformation.CANGLE[1],
					prediction_transformation.CANGLE[2],
					prediction_transformation.TVECT[0], prediction_transformation.TVECT[1],
					prediction_transformation.TVECT[2],
					predictions.N, predictions.SPACING);
		}

		/*
		 * Returns the ZDOCK score assigned to the prediction number specified by the index
		 */
		double get_score(size_t transformation_number)
		{
			// get the specific transformation...
			zdock prediction_transformation = predictions.ZD[transformation_number];
			// and return the score
			return prediction_transformation.SCORE;
		}

		/*
		 * Return the size of the zdata instance stored
		 */
		size_t get_size()
		{
			return predictions.ZD.size();
		}

		/*
		 * Prints all the information corresponding to each transformation
		 * to the output stream provided
		 */
		void print(ostream& output_stream)
		{
			// print the general information first
			output_stream << predictions.N << "\t" << predictions.SPACING << endl;
			output_stream << predictions.RANDVECT[0] << "\t" << predictions.RANDVECT[1] << "\t" << predictions.RANDVECT[2] << endl;
			output_stream << predictions.RVECT[0] << "\t" << predictions.RVECT[1] << "\t" << predictions.RVECT[2] << endl;
			output_stream << predictions.LVECT[0] << "\t" << predictions.LVECT[1] << "\t" << predictions.LVECT[2] << endl;

			// print each transformation
			for(size_t current_index = 0; current_index < predictions.ZD.size(); current_index++)
			{
				zdock current = predictions.ZD[current_index];
				output_stream << current.CANGLE[0] << "\t" << current.CANGLE[1] << "\t" << current.CANGLE[2] << "\t"
						<< current.TVECT[0] << "\t" << current.TVECT[1] << "\t" << current.TVECT[2] << "\t"
						<< current.SCORE << endl;
			}
		}

		/*
		 * Inverting a ZDOCK transformation consists of changing the sign of the rotation
		 * angle and the x,y,z translation values. We apply the transformation for each
		 * of the predictions stored. Additionally, the starting positions of the ligand
		 * and receptor must be switched, because the roles are reversed
		 * Note: The initial rotation of the ligand (set with -S parameter in ZDOCK) is not used
		 * to generate decoys, thus it's not taken into account
		 */
		void invert()
		{
			/* create a copy of the original receptor translation */
			double oldreceptorx = predictions.RVECT[0];
			double oldreceptory = predictions.RVECT[1];
			double oldreceptorz = predictions.RVECT[2];
			/* and now overwrite the receptor initial x,y,z */
			predictions.RVECT[0] = predictions.LVECT[0];
			predictions.RVECT[1] = predictions.LVECT[1];
			predictions.RVECT[2] = predictions.LVECT[2];
			/* and now set the original ligand x,y,z to the old receptor values */
			predictions.LVECT[0] = oldreceptorx;
			predictions.LVECT[1] = oldreceptory;
			predictions.LVECT[2] = oldreceptorz;
	
			/* Invert all of the transformations */
			for(size_t current_index = 0; current_index < predictions.ZD.size(); current_index++)
			{
				zdock current = predictions.ZD[current_index];
				
				/* change the sign of the rotation and translation parameters */
				predictions.ZD[current_index].CANGLE[0] *= -1;
				predictions.ZD[current_index].CANGLE[1] *= -1;
				predictions.ZD[current_index].CANGLE[2] *= -1;

				predictions.ZD[current_index].TVECT[0] *= -1;
				predictions.ZD[current_index].TVECT[1] *= -1;
				predictions.ZD[current_index].TVECT[2] *= -1;
			}
		
		}
	};
#endif
