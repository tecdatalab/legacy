#ifndef _MRF_SAMPLER_H_
#define _MRF_SAMPLER_H_

#include <map>
#include <utility>
#include <vector>

#include "../em/density_lattice.h"

#include "mrf_features.h"
#include "pdb.h"
#include "transformable_pdb.h"

using std::pair;
using std::vector;

// Instances of this class are used to generate data files that sample different
// positions of PDB units fittted into an EM map. These files will contain
// features that are later used to train models using an MRF framework.
// The sampling process starts by finding the best position for each unit,
// with the help of the c-alpha trace, and then samples a small number of
// positions close to the best location (short_rotation_* parameters) as well
// as a more widespread sampling (using long_rotation_*).
// For each unit X, a file called X.mrf is created with CC, RMSD and Overlap
// features, and for each neighbor pair X-Y an analog X-Y.mrf is created
// with pairwise shape and physics-based features.
// The small sampling is intended to generate close to native data points,
// while the long sampling should generate a varied enough dataset by sampling
// across the whole rotational space.
class mrf_sampler {
 public:
  // pdb_unit_files: Each atomic description of the units e.g. A.pdb, B.pdb, etc
  // calpha_trace_file: Best arrangement of the PDB units that can be fitted
  //   into the complete_map_file. All the pdb_unit_files can be superimposed
  //   onto this description of the c-alpha trace
  // em_unit_files: Analog to pdb_unit_files, but MRC EM maps instead
  // complete_map_file: Map that matches the c-alpha trace
  // contour_level: Contour value applied to all EM maps, normally from EMDB,
  // unit_labels: Used to later identify the neighbor relations. Must be in the
  //   same order as the PDB and EM units provided
  // neighbors: edges between units, used to generate pairwise features
  // short_rotation_degree_range: from the initial structure +- these number
  //   of degrees
  // short_rotation_step: starting from -(short_rotation_degree_range) add this
  //   step value until getting to +(short_rotation_degree_range)
  // short_translation_range/step: analog to the rotation range, but this is
  //   +- angstroms in translations following the step size parameter
  // long_rotation_step: same thing but without the range parameter because we
  //   sample all 360 degrees
  // long_translation_range/step: same thing but for the long type of sampling.
  // output_prefix: adds a prefix to sample files
  // pdb_output_prefix: if non-empty triggers the generation of PDB files
  //   (with filenames having this prefix) for each singleton transformation.
  //   If provided, pairwise files are not generated.
  mrf_sampler(vector<string> pdb_unit_files, string calpha_trace_file,
      vector<string> em_unit_files, string complete_map_file,
      double contour_level, vector<string> unit_labels,
      vector<pair<string, string> > neighbors,
      double short_rotation_degree_range, double short_rotation_step,
      double short_translation_range, double short_translation_step,
      double long_rotation_step, double long_translation_range,
      double long_translation_step, string output_prefix,
      string pdb_output_prefix);
  // Density lattices are allocated dynamically, release them.
  ~mrf_sampler();
  // Main function that starts the process of sampling and generating the
  // one and two-body feature files.
  // swapping_enabled: true when alternate conformations for each position
  // need to be generated. e.g. if there are 3 units A,B,C with swapping_enabled
  // conformations of B and C placed at A's center are generated, but no
  // conformations with A in its intended position are created.
  void generate_feature_files(bool swapping_enabled);
 private:
  // Loads all em map units and the reference map into instance variables.
  void initialize_em_maps(string reference_map_filename,
                          vector<string> em_unit_filenames);
  // Set the aligned_units instance member by transforming the input PDBs to
  // optimally align to the c-alpha trace.
  void initialize_pdb_units(vector<string> pdb_filenames,
                            vector<string> labels,
                            pdb& calpha_template);
  // Called at the end of initialize_pdb_units in order to create all the
  // point_transformations needed in case unit swapping is enabled during
  // the sampling process.
  void initialize_swap_transformations();
  // These two methods generate the CC, overlap and RMSD values.
  // The first one iterates over all units and transformation and calls
  // the second one that is in charge of computing the values for a single
  // unit and transformation.
  // swapping_enabled: as described in generate_feature_files
  void generate_one_body_files(double alpha_center, double beta_center,
      double gamma_center, double alpha_plus_minus, double beta_plus_minus,
      double gamma_plus_minus, double translation_plus_minus,
      double rotation_step, double translation_step, bool append_to_file,
      bool swapping_enabled);
  // base_transformation: The x,y,z,alpha,beta,gamma sampled
  // applied_transformation: Contains base_transformation but it also
  // has additional adjustments like moves to the centroid
  mrf_one_body_features calculate_one_body_features(size_t unit_index,
      point_transformation& base_transformation,
      point_transformation_sequence& applied_transformation);
  // Generate a single PDB file by modifying protein[unit_index] with
  // the transformation given and storing it in
  // <pdb_prefix>-<unit_label>-<file_sequence_number>.pdb
  void generate_singleton_pdb(size_t unit_index,
      point_transformation_sequence& applied_transformation, string pdb_prefix,
      size_t file_sequence_number);
  // Iteration over all edges in the mrf and 6D transformation for each
  // of the 2 units involved. Similar to the previous methods, the second
  // one computes the features for a particular transformation of the 2
  // vertices in the graph
  // swapping_enabled: as described in generate_feature_files
  void generate_two_body_files(double alpha_center, double beta_center,
      double gamma_center, double alpha_plus_minus, double beta_plus_minus,
      double gamma_plus_minus, double translation_plus_minus,
      double rotation_step, double translation_step, bool append_to_file,
      bool swapping_enabled);
  // left_pdb and right_pdb have been already transformed. The
  // point_transformations are provided as parameters to correctly
  // initialize the mrf_two_body_features instance.
  // Note: base vs applied transformations are the analog to what
  // is done with the one-body features
  mrf_two_body_features calculate_two_body_features(
      size_t left_index, size_t right_index,
      point_transformation& base_left_transformation,
      point_transformation& base_right_transformation,
      point_transformation_sequence& applied_left_transformation,
      point_transformation_sequence& applied_right_transformation);
  // Converts unit indices to actual hash table keys to retrieve the appropriate
  // transformation that moves move_from_unit to move_to_unit
  point_transformation get_swap_transformation(size_t move_from_unit,
      size_t move_to_unit);
  // If no swapping is done, a single base_transformation is applied to each unit
  // with adjustments done to move to the centroid first.
  // If swapping is enabled then the function returns not just a single
  // transformation_sequence but (total_units - 1), each representing the swap of
  // one unit by another.
  // Return value: If no swapping, a vector with one pair, where the unit_index
  // is the same as the one input. If swapping is enabled, unit_indices will be
  // all except the unit_index input along with the transformation.
  vector<pair<size_t, point_transformation_sequence> >
      get_all_adjusted_transformations(size_t unit_index,
                                       point_transformation base_transformation,
                                       bool swapping_enabled);
  // Since transformations are calculated from each unit's centroid,
  // this method adds three point_transformations to the sequence:
  // a translation to the centroid, the actual transformation and 
  // a translation that negates the original centroid move.
  void add_centroid_adjusted_transformation(size_t unit_index,
      point_transformation transformation,
      point_transformation_sequence* sequence);
  // It should be a PDB ID referenced from the EMDB entry that provides
  // the best fitted superimposition of the c-alpha atoms into the map
  pdb calpha_trace;
  // Matches calpha_trace, loaded from an EM map obtained directly from
  // the EMDB.
  density_lattice* reference_em;
  // Normally the same value for the reference map found in the EMDB entry.
  double contour_level;
  // Full atomic PDB fragments that have been transformed to align as well
  // as possible to the c-alpha trace. These are considered to be the starting
  // atomic positions for sampling purposes.
  vector<transformable_pdb> aligned_units;
  // Transformations applied to input PDBs to optimally align them to
  // the calpha_trace
  vector<point_transformation> alignment_transformations;
  // Each density lattice represents an EM map, in the same order as the units.
  // These maps will normally be created from simulated EM map files.
  vector<density_lattice*> em_units;
  // For each map and the contour_level provided, we precompute the boundaries
  // of the EM surfaces
  vector<density_vector> em_min_boundaries;
  vector<density_vector> em_max_boundaries;
  // MRF Graph description. Node names and edges.
  vector<string> unit_labels;
  vector<pair<string,string> > neighbors;
  // Keep track of sequence numbers used to generate PDB files. This is just
  // to make the output files have the same number in their filename as the
  // row number each transformation has in the sampling output files.
  vector<int> transformation_counts;
  // Quick reference to determine what unit index corresponds to each label
  map<string, size_t> label_indices;
  // Keys are X-Y where X and Y are unit labels. Values associated to keys
  // are translations that move X's centroid to Y's centroid.
  // This used when generating samples for units that are swapped, becoming
  // negative examples.
  map<string, point_transformation> swap_transformations;
  // Sampling parameters.
  double short_rotation_degree_range;
  double short_rotation_step;
  double short_translation_range;
  double short_translation_step;
  double long_rotation_step;
  double long_translation_range;
  double long_translation_step;
  // Files created have this prefix.
  string output_prefix;
  // Same for PDB files
  string pdb_output_prefix;
};
#endif
