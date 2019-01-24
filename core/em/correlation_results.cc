#include "correlation_results.h"

#include "../proteinstructure/point_transformation.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;

correlation_result::correlation_result() : is_center_set(false) {
}

void correlation_result::calculate_transformed_center(
    double center_x, double center_y, double center_z) {
  point_transformation transformation(alpha_z, beta_y, gamma_x, x, y, z);
  transformation.transform(center_x, center_y, center_z,
      &transformed_center_x, &transformed_center_y, &transformed_center_z);
  this->is_center_set = true;
}

double* correlation_result::coordinates()
{
  double* joint_coordinates = new double[3];
  joint_coordinates[0] = transformed_center_x;
  joint_coordinates[1] = transformed_center_y;
  joint_coordinates[2] = transformed_center_z;

  return joint_coordinates;
}

bool correlation_results::compare(correlation_result result1,
								correlation_result result2)
{
	return result1.correlation > result2.correlation ||
			// although it's not common, the tiebraker is the overlap
			(result1.correlation == result2.correlation &&
			result1.overlap > result2.overlap);
}

void correlation_results::read(std::string filename)
{
	// These are just a tmp variables to read each result
	// before pushing into the collection
  char results_line[1024];
	ifstream results_file(filename.c_str());
  if(!results_file.is_open()) {
    // just a generic exception number
    cerr << "Could not open results file" << endl;
    throw 1;
  }
	while(results_file.good())
	{
	  correlation_result current_result;
    results_file.getline(results_line, 1024);
    istringstream line_tokenizer(results_line);
    // We expect at least x,y,z, three angles, cc and overlap.
		line_tokenizer >> current_result.x >> current_result.y >> current_result.z
					>> current_result.gamma_x >> current_result.beta_y
					>> current_result.alpha_z >> current_result.correlation
					>> current_result.overlap;
    // If there's still elements left, we assume it will be the transformed_center
    if(line_tokenizer.good()) {
      line_tokenizer >> current_result.transformed_center_x
          >> current_result.transformed_center_y
          >> current_result.transformed_center_z;
      current_result.is_center_set = true;
    }
		push_back(current_result);
	}
	results_file.close();
}

void correlation_results::write(std::string filename)
{
	ofstream results_file(filename.c_str());
    if(!results_file.is_open()) {
        // just a generic exception number
        cerr << "Could not create results file" << endl;
        throw 1;
    }

	// Before writing, sort the results
	sort(begin(), end(), compare);	

	correlation_results::iterator results_it;
	for(results_it = begin(); results_it < end(); results_it++)
	{
		correlation_result current_result = *results_it;
    results_file.unsetf(std::ios_base::floatfield);
		results_file << current_result.x << " " << current_result.y << " "
					<< current_result.z << " " << current_result.gamma_x << " "
					<< current_result.beta_y << " "
					<< current_result.alpha_z << " ";
		results_file << std::fixed << current_result.correlation << " "
					       << current_result.overlap;
    if(current_result.is_center_set) {
      results_file << " " << std::fixed << current_result.transformed_center_x
                   << " " << current_result.transformed_center_y
                   << " " << current_result.transformed_center_z;
    }
    results_file << endl;
	}

	results_file.close();
}
