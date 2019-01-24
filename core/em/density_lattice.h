#ifndef _DENSITY_LATTICE_H_
#define _DENSITY_LATTICE_H_

#include "../proteinstructure/point_transformation.h"

#include "correlation_results.h"
#include "point_lattice.h"

// This is in voxel space, namely the coordinates go from 0 to
// (numberOfVoxels - 1).
typedef struct
{
	unsigned long x, y, z;
	float density;
} density_voxel;

// This is the actual x,y,z coordinate in real space for a density
// vector. To get the x value for this vector we just multiply
// the x in a density voxel by the actual voxel_length_x
typedef struct
{
  double x, y, z;
  float density;
} density_vector;

/*
 * Simplified representation of a density map
 * All values are assumed to be float and
 * the description of the number of voxels and
 * size of each voxel is stored as part of this instance
 */

class density_lattice : public point_lattice
{
 private:
  // The min and max boundaries only give two points out of the 8 that define
  // the bounding box that contains the density surface.
  // This method simply infers the other points from the coordinates
  // in min and max and puts all 8 in the vector returned.
  //
  // This method is used as an auxiliary by get_boundaries_after_transformation
  vector<density_vector> calculate_all_corners(const density_vector& min,
                                               const density_vector& max);
 public:
    density_lattice(unsigned long in_nx, unsigned long in_ny,
                    unsigned long in_nz, float in_xorigin, float in_yorigin,
                    float in_zorigin, float in_voxel_xlen, float in_voxel_ylen,
                    float in_voxel_zlen) :
        point_lattice(in_nx, in_ny, in_nz, in_xorigin, in_yorigin, in_zorigin,
            in_voxel_xlen, in_voxel_ylen, in_voxel_zlen) {}
    // Scans all voxels and returns min and max values
    void get_min_max_densities(float* min, float* max);
	  // Iterates over the whole density lattice and returns (by reference)
	  // The values of the lowest x,y,z coordinates that have a density
	  // >= to density threshold, and conversely for the highest x,y,z.
    // The coordinates are in real space (not voxel space).
	  void get_boundaries_for_threshold(float density_threshold,
		                  								density_vector* min_bound,
										                  density_vector* max_bound);
    // Analog to the previous one but the boundaries are given in voxel space.
	  void get_voxel_boundaries_for_threshold(float density_threshold,
                                            density_voxel* min_bound,
                                            density_voxel* max_bound);
    // The previous get boundary methods assume that the map is in its
    // original orientation, however, once a transformation is applied,
    // we need to recalculate what the min and max values along x, y and z
    // axis are. This is because we want to compare the least amount of
    // density voxels possible when compute the correlation and overlap
    // between two maps.
    // 
    // min_before and max_before should be obtained by calling
    // get_boundaries_for_threshold first. They are passed as parameters instead
    // of computing them within the method to avoid recalculations, since
    // getting the boundaries implies going through all the voxels.
    void get_boundaries_after_transformation(
        point_transformation_sequence transform,
        density_vector min_before, density_vector max_before,
        density_vector* min_after, density_vector* max_after);
    // Assuming both lattices are comparable (e.g. same voxel size)
    // This method will match, for all valid [x,y,z],
    // this->densities[x+offset_x,y+offset_y,z+offset_z] (a) with
    // compared_lattice[x,y,z] (b)
    // and computes the correlation: sum(a*b)/sqrt(sum(a^2)*sum(b^2))
    double calculate_correlation(density_lattice* compared_lattice,
                                 unsigned long offset_x,
                                 unsigned long offset_y,
                                 unsigned long offset_z,
                								 double alpha_z, double beta_y, double gamma_x,
                                 float density_threshold,
                                 unsigned long* overlap);
  // The previous version is legacy code. This new version assumes that
  // real space coordinates are provided as bounds to sample this lattice.
  // This methods samples from min_bound to max bound (with a step size of 
  // "voxel length").
  // For each coordinate it will use the value at x,y,z in this lattice
  // and the inverse_transform(x,y,z) to sample the compared_lattice
  // (Note that applying the inverse to find individual points
  // is the same as trying to transform the whole lattice).
  //
  // The correlation is computed as: sum(a*b)/sqrt(sum(a^2)*sum(b^2))
  // where a and b correspond to this_lattice and compared_lattice
  // 
  // The overlap is: numerator->voxels where both densities >= threshold
  //                 denominator->voxels in compared_lattice >= threshold
  double calculate_correlation(density_lattice& compared_lattice,
      point_transformation_sequence& compared_transformation,
      density_vector min_bound, density_vector max_bound,
      float density_threshold, double* overlap);
	// This method calls calculate correlation for different values
	// of offset_x/y/z. It will start from 0,0,0 all the way up to
	// max_x/y/z, adding step_size to each coordinate iteratively
	// For each translation, it will also try rotations from 0 up to
	// 360 degrees at most, using rotation_step for each increment
	correlation_results get_significant_correlations(
		density_lattice* compared_lattice, int step_size, double rotation_step,
		float density_threshold, float correlation_threshold,
		float overlap_threshold);
	// Returns the number of voxels (times the voxel volume)
	// that have a density >= than the value provided
	double calculate_volume(float density_threshold);
	// Returns the number of voxels that have a density >= the threshold
	// (Not multiplied by the voxel volume)
	unsigned long count_voxels_over_threshold(float density_threshold);
};

#endif
