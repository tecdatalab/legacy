#ifndef _CORRELATION_RESULTS_H_
#define _CORRELATION_RESULTS_H_

#include <string>
#include <vector>

// Correlation results are just the translation/rotation
// transformation applied to a density map and the
// subsequent correlation coefficient and overlap with respect
// to some other map
class correlation_result
{
 public:
	double x, y, z;
	double alpha_z, beta_y, gamma_x;
	float overlap;
	double correlation;
  // These represent the center coordinates of an EM volume, for a
  // particular isosurface threshold. They are not part of the standard
  // result output to files, but they are used to cluster the results,
  // since the clustering should compute distances from the volume centers.
  double transformed_center_x, transformed_center_y, transformed_center_z;
  // Upon creation this flag is set to false. After the first call to
  // calculate_transformed_center it's set to true and kept that way
  // for the remainder of the instance's life.
  bool is_center_set;
  correlation_result();
  // This method receives the volume center, before applying the
  // transformation stored in this instance. It will take the coordinates,
  // and apply the transformation in order to store them in the
  // transformed_center_? instance variables.
  void calculate_transformed_center(double center_x, double center_y,
                                    double center_z);
  // Used by kmeans to generically get the x,y,z values to
  // compute distances and centroids. In this case they represent
  // the volume center coordinates. It is a requirement for
  // calculate_transformed_center to have been called before.
  double* coordinates();
};

/*
 * This is a collection of correlation_result.
 * It is in itself a vector but it has methods to write or read
 * from a file
 */
class correlation_results : public std::vector<correlation_result>
{
  private:
	// Before writing files, results are sorted by correlation
	// and overlap so we need a comparison function to call
	// sort on
	static bool compare(correlation_result result1, correlation_result result2);
  public:
	// No constructor necessary, just the default initialization
	// performed automatically

	// Read from a file. Assume the following space-separated values
	// are provided in each line (one line per correlation_result):
	// translation_x translation_y translation_z
	// rotation_x rotation_y rotation_z
	// correlation overlap
	void read(std::string filename);
	// Writes a header with column names and then writes one line per
	// correlation_result, as described in the read method
	void write(std::string filename);	
};
#endif
