#include "mrf_sampler.h"

#include <fstream>
#include <utility>

#include "../em/mrc.h"

#include "align.h"
#include "atom.h"
#include "mrf_features.h"
#include "point_transformation.h"
#include "rmsd.h"
#include "scoring.h"
// Included only because of the soroban_score predefined weights
#include "score_weights_options.h"
#include "soroban_score.h"

mrf_sampler::mrf_sampler(vector<string> pdb_unit_files,
    string calpha_trace_file, vector<string> em_unit_files,
    string complete_map_file, double contour_level, vector<string> unit_labels,
    vector<pair<string, string> > neighbors,
    double short_rotation_degree_range, double short_rotation_step,
    double short_translation_range, double short_translation_step,
    double long_rotation_step, double long_translation_range,
    double long_translation_step, string output_prefix,
    string pdb_output_prefix) :
        contour_level(contour_level),
        unit_labels(unit_labels),
        neighbors(neighbors),
        short_rotation_degree_range(short_rotation_degree_range),
        short_rotation_step(short_rotation_step),
        short_translation_range(short_translation_range),
        short_translation_step(short_translation_step),
        long_rotation_step(long_rotation_step),
        long_translation_range(long_translation_range),
        long_translation_step(long_translation_step),
        output_prefix(output_prefix), pdb_output_prefix(pdb_output_prefix) {
  // 1. Initialize the calpha trace, used to align and later to compute RMSDs.
  read_protein(calpha_trace_file, calpha_trace);
  // 2. Load the  EM maps.
  initialize_em_maps(complete_map_file, em_unit_files);
  // 3. Set the initial position for the PDB files (initial optimal placement).
  initialize_pdb_units(pdb_unit_files, unit_labels, calpha_trace);
  // 4. For later use, initialize the mapping between labels and the unit
  // index they represent
  for (size_t unit_index = 0; unit_index < unit_labels.size(); unit_index++) {
    string label = unit_labels[unit_index];
    label_indices[label] = unit_index;
  }
  // Intiialize all transformation sequence numbers to zero
  transformation_counts.insert(transformation_counts.begin(),
                               pdb_unit_files.size(), 0);
}

mrf_sampler::~mrf_sampler() {
  if (reference_em) {
    delete reference_em;
  }
  for (size_t map_index = 0; map_index < em_units.size(); map_index++) {
    density_lattice* map = em_units[map_index];
    if (map) {
      delete map;
    }
  }
}

void mrf_sampler::generate_feature_files(bool swapping_enabled) {
  // The short sampling creates a new file with the short_xxx parameters and
  // then a second call appends using the long sampling parameters
  generate_one_body_files(0, 0, 0, short_rotation_degree_range,
      short_rotation_degree_range, short_rotation_degree_range,
      short_translation_range, short_rotation_step, short_translation_step,
      false, swapping_enabled);
  generate_one_body_files(179, -1, 179, 179, 89, 179, long_translation_range,
      long_rotation_step, long_translation_step, true, swapping_enabled);

  // Ignore two body files if PDB's are expected. This is just to save time if
  // only the PDB's are needed.
  if (pdb_output_prefix.empty()) {
    generate_two_body_files(0, 0, 0, short_rotation_degree_range,
        short_rotation_degree_range, short_rotation_degree_range,
        short_translation_range, short_rotation_step, short_translation_step,
        false, swapping_enabled);
    generate_two_body_files(179, -1, 179, 179, 89, 179, long_translation_range,
        long_rotation_step, long_translation_step, true, swapping_enabled);
  }
}

void mrf_sampler::initialize_em_maps(string reference_map_filename,
    vector<string> em_unit_filenames) {
  mrc_reader reference_reader(reference_map_filename);
  reference_em = reference_reader.get_densities();
  
  for (size_t file_index = 0; file_index < em_unit_filenames.size();
       file_index++) {
    mrc_reader map_reader(em_unit_filenames[file_index]);
    em_units.push_back(map_reader.get_densities());
  } 
}

void mrf_sampler::initialize_pdb_units(vector<string> pdb_filenames,
                                       vector<string> labels,
                                       pdb& calpha_template) {
  for(size_t unit_index = 0; unit_index < pdb_filenames.size(); unit_index++) {
    pdb pdb_unit;
    read_protein(pdb_filenames[unit_index], pdb_unit);
    // Create the optimally aligned representation
    vector<atom> input_matches;
    vector<atom> trace_matches;

    pdb_unit.get_matching_calphas(&calpha_template, true, // ignores chain ID
                                  &input_matches, &trace_matches);
    if(input_matches.size() != trace_matches.size() ||
       (!input_matches.size())) {
      cerr << "Error aligning unit number " << (unit_index+1)
           << ". All c-alpha atoms in the input units must match at least one "
           << "atom in the c-alpha trace." << endl;
      exit(EXIT_FAILURE);
    }
    vector<atom> aligned_output;
    point_transformation alignment_transformation;
    double rmsd = align_all_and_transform(
        &input_matches, &trace_matches, &pdb_unit.atoms, &aligned_output,
        &alignment_transformation);
    cout << "Unit " << labels[unit_index] << " aligned with an RMSD of "
         << rmsd << " angstroms" << endl;
    transformable_pdb aligned_unit(labels[unit_index], aligned_output);
    aligned_unit.precompute_centroid();
    double centroid_x, centroid_y, centroid_z;
    aligned_unit.get_centroid(&centroid_x, &centroid_y, &centroid_z);
    cout << "Aligned centroid: " << centroid_x << " " << centroid_y << " "
         << centroid_z << endl;

    aligned_units.push_back(aligned_unit);
    alignment_transformations.push_back(alignment_transformation);
    // Precompute the boundaries for the particular threshold used by this
    // instance of mrf_sampler.
    density_vector min, max;
    density_lattice* current = em_units[unit_index];
    current->get_boundaries_for_threshold(contour_level, &min, &max);
    em_min_boundaries.push_back(min);
    em_max_boundaries.push_back(max);

    // Testing the CC and OVR between the best alignment and the map
    point_transformation_sequence seq;
    seq.add(alignment_transformation);

    density_vector tmin, tmax;
    current->get_boundaries_after_transformation(seq, min, max, &tmin, &tmax);

    cout << "Initial EM boundaries: Min " << min.x << " " << min.y << " "
         << min.z << " Max " << max.x << " " << max.y << " " << max.z << endl;
    cout << "Optimal EM boundaries: Min " << tmin.x << " " << tmin.y << " "
         << tmin.z << " Max " << tmax.x << " " << tmax.y << " " << tmax.z
         << endl;

    double test_cc, test_overlap;
    test_cc = reference_em->calculate_correlation(*current, seq, tmin, tmax,
                                             contour_level, &test_overlap);
    cout << "CC " << test_cc << " Overlap " << test_overlap << "\n";

    initialize_swap_transformations();
  }
}

void mrf_sampler::initialize_swap_transformations() {
  for (size_t move_from_index = 0; move_from_index < aligned_units.size();
       move_from_index++) {
    for (size_t move_to_index = move_from_index + 1;
         move_to_index < aligned_units.size(); move_to_index++) {
      // Get both centroids
      double from_centroid_x, from_centroid_y, from_centroid_z;
      double to_centroid_x, to_centroid_y, to_centroid_z;
      aligned_units[move_from_index].get_centroid(&from_centroid_x,
          &from_centroid_y, &from_centroid_z);
      aligned_units[move_to_index].get_centroid(&to_centroid_x,
          &to_centroid_y, &to_centroid_z);
      // We generate both translations to move either unit to the
      // other's centroid
      string from_label = unit_labels[move_from_index];
      string to_label = unit_labels[move_to_index];
      string from_to_swap_key = from_label + "-" + to_label;
      string to_from_swap_key = to_label + "-" + from_label;
      double translation_x = -(from_centroid_x - to_centroid_x);
      double translation_y = -(from_centroid_y - to_centroid_y);
      double translation_z = -(from_centroid_z - to_centroid_z);
      // The zeroes represent 0 degree rotations
      swap_transformations[from_to_swap_key] = point_transformation(0, 0, 0,
          translation_x, translation_y, translation_z);
      swap_transformations[to_from_swap_key] = point_transformation(0, 0, 0,
          -translation_x, -translation_y, -translation_z);
    }
  }
}

void mrf_sampler::generate_one_body_files(double alpha_center,
    double beta_center, double gamma_center, double alpha_plus_minus,
    double beta_plus_minus, double gamma_plus_minus,
    double translation_plus_minus, double rotation_step,
    double translation_step, bool append_to_file, bool swapping_enabled) {
  for (size_t unit_index = 0; unit_index < aligned_units.size(); unit_index++) {
    point_transformation_generator transformations(alpha_center, beta_center,
        gamma_center, alpha_plus_minus, beta_plus_minus, gamma_plus_minus,
        0, 0, 0, // translation centers
        translation_plus_minus, translation_plus_minus, translation_plus_minus,
        rotation_step, translation_step);
    // This method is called twice per file, the first time will truncate it,
    // the second will append.
    ios_base::openmode mode = append_to_file ? std::ofstream::app :
                                               std::ofstream::trunc;
    ofstream sampling_results_file(
        (output_prefix + "-" + unit_labels[unit_index] + ".mrf").c_str(), mode);
    if (!append_to_file) {
      mrf_one_body_features::write_header(sampling_results_file);
    }
    while (!transformations.is_finished()) {
      point_transformation base_transformation = transformations.next();
      vector<pair<size_t, point_transformation_sequence> > samples =
          get_all_adjusted_transformations(unit_index, base_transformation,
                                           swapping_enabled);
      for (size_t sample_index = 0; sample_index < samples.size();
           sample_index++) {
        size_t applied_index = samples[sample_index].first;
        point_transformation_sequence applied_transformation =
            samples[sample_index].second;

        mrf_one_body_features features = calculate_one_body_features(
            applied_index, base_transformation, applied_transformation);
        features.write(sampling_results_file);
        if (!pdb_output_prefix.empty()) {
          generate_singleton_pdb(applied_index, applied_transformation,
              pdb_output_prefix, transformation_counts[unit_index]);
        }
        transformation_counts[unit_index]++;
      }
    }
    sampling_results_file.close();
  }
}

mrf_one_body_features mrf_sampler::calculate_one_body_features(
    size_t unit_index, point_transformation& base_transformation,
    point_transformation_sequence& applied_transformation) {
  transformable_pdb current_unit = aligned_units[unit_index];
  current_unit.apply_point_transformation_sequence(&applied_transformation);

  vector<atom> unit_matches;
  vector<atom> trace_matches;
  // PDB related computations
  current_unit.get_matching_calphas(&calpha_trace, true,
      &unit_matches, &trace_matches);
  double rmsd = calculate_allatom_rmsd(unit_matches, trace_matches);

  // EM map related computations
  density_vector tmin, tmax;
  density_vector min = em_min_boundaries[unit_index];
  density_vector max = em_max_boundaries[unit_index];
  density_lattice* current = em_units[unit_index];
  point_transformation_sequence seq;
  // Since the EM map hasn't been pre-aligned we first need to move it to
  // the alignment site
  seq.add(alignment_transformations[unit_index]);
  seq.add(applied_transformation);
  current->get_boundaries_after_transformation(seq, min, max, &tmin, &tmax);

  double overlap;
  double cc = reference_em->calculate_correlation(*current, seq, tmin, tmax,
      contour_level, &overlap);

  mrf_one_body_features features(base_transformation.get_alpha(),
      base_transformation.get_beta(), base_transformation.get_gamma(),
      base_transformation.get_tx(), base_transformation.get_ty(),
      base_transformation.get_tz(), rmsd, cc, overlap, unit_labels[unit_index]);
  
  return features;
}

void mrf_sampler::generate_singleton_pdb(size_t unit_index,
    point_transformation_sequence& applied_transformation, string pdb_prefix,
    size_t file_sequence_number) {
  transformable_pdb current_unit = aligned_units[unit_index];
  current_unit.apply_point_transformation_sequence(&applied_transformation);
  write_complex(pdb_prefix + "-" + unit_labels[unit_index],
                file_sequence_number, current_unit.atoms);
}

void mrf_sampler::generate_two_body_files(double alpha_center,
    double beta_center, double gamma_center, double alpha_plus_minus,
    double beta_plus_minus, double gamma_plus_minus,
    double translation_plus_minus, double rotation_step,
    double translation_step, bool append_to_file, bool swapping_enabled) {
  for (size_t edge_index = 0; edge_index < neighbors.size(); edge_index++) {
    pair<string,string> edge = neighbors[edge_index];
    size_t left_index = label_indices[edge.first];
    size_t right_index = label_indices[edge.second];

    point_transformation_generator left_transformations(alpha_center,
        beta_center, gamma_center, alpha_plus_minus, beta_plus_minus,
        gamma_plus_minus, 0, 0, 0, // translation centers
        translation_plus_minus, translation_plus_minus, translation_plus_minus,
        rotation_step, translation_step);
    // This method is called twice per file, the first time will truncate it,
    // the second will append.
    ios_base::openmode mode = append_to_file ? std::ofstream::app :
                                               std::ofstream::trunc;
    ofstream sampling_results_file(
        (output_prefix + "-" + unit_labels[left_index] + "-"
                       + unit_labels[right_index] + ".mrf").c_str(), mode);
    if (!append_to_file) {
      mrf_two_body_features::write_header(sampling_results_file);
    }
    while (!left_transformations.is_finished()) {
      point_transformation base_left_transformation =
          left_transformations.next();
      vector<pair<size_t, point_transformation_sequence> > left_samples =
          get_all_adjusted_transformations(left_index, base_left_transformation,
                                           swapping_enabled);
      for (size_t left_sample_index = 0; left_sample_index < left_samples.size();
           left_sample_index++) {
        size_t applied_left_index = left_samples[left_sample_index].first;
        point_transformation_sequence applied_left_transformation =
            left_samples[left_sample_index].second;
        // Iterate over the right PDB transformations
        point_transformation_generator right_transformations(alpha_center,
            beta_center, gamma_center, alpha_plus_minus, beta_plus_minus,
            gamma_plus_minus, 0, 0, 0, // translation centers
            translation_plus_minus, translation_plus_minus,
            translation_plus_minus, rotation_step, translation_step);
        while (!right_transformations.is_finished()) {
          point_transformation base_right_transformation =
              right_transformations.next();
          vector<pair<size_t, point_transformation_sequence> > right_samples =
              get_all_adjusted_transformations(right_index,
                  base_right_transformation, swapping_enabled);
          for (size_t right_sample_index = 0;
               right_sample_index < right_samples.size();
               right_sample_index++) {
            size_t applied_right_index =
                right_samples[right_sample_index].first;
            // When swapping is enabled we want to prevent the same unit to be
            // located in two places
            if (applied_left_index != applied_right_index) {
              point_transformation_sequence applied_right_transformation = 
                  right_samples[right_sample_index].second;
              mrf_two_body_features features = calculate_two_body_features(
                  applied_left_index, applied_right_index,
                  base_left_transformation, base_right_transformation,
                  applied_left_transformation, applied_right_transformation);
              features.write(sampling_results_file);
            }
          }
        }
      }
    }
    sampling_results_file.close();
  }
}

mrf_two_body_features mrf_sampler::calculate_two_body_features(
    size_t left_index, size_t right_index,
    point_transformation& base_left_transformation,
    point_transformation& base_right_transformation,
    point_transformation_sequence& applied_left_transformation,
    point_transformation_sequence& applied_right_transformation) {
  transformable_pdb left_pdb = aligned_units[left_index];
  transformable_pdb right_pdb = aligned_units[right_index];
  // Apply translations wrt the centroids
  left_pdb.apply_point_transformation_sequence(&applied_left_transformation);
  right_pdb.apply_point_transformation_sequence(&applied_right_transformation);

  // To calculate the RMSD merge both sets of atoms to align them to
  // the c-alpha trace
  vector<atom> transformed_matches;
  vector<atom> trace_matches;
  pdb* both_pdb = pdb::merge(left_pdb, right_pdb);
  both_pdb->get_matching_calphas(&calpha_trace, true,
      &transformed_matches, &trace_matches);

  double rmsd = calculate_allatom_rmsd(transformed_matches, trace_matches);

  // Calculate the physics score terms (it requires two separate vectors)
  vector<vector<atom> > decoy;
  decoy.push_back(left_pdb.atoms);
  decoy.push_back(right_pdb.atoms);

  // Zero means that we're using the weights that use all 12 terms
  soroban_score score = compute_energy(decoy, TemplateWeights[0]);
  mrf_two_body_features features(base_left_transformation.get_alpha(),
      base_left_transformation.get_beta(), base_left_transformation.get_gamma(),
      base_left_transformation.get_tx(), base_left_transformation.get_ty(),
      base_left_transformation.get_tz(), base_right_transformation.get_alpha(),
      base_right_transformation.get_beta(),
      base_right_transformation.get_gamma(),
      base_right_transformation.get_tx(), base_right_transformation.get_ty(),
      base_right_transformation.get_tz(), rmsd, score, unit_labels[left_index],
      unit_labels[right_index]);
  delete both_pdb;
  return features;
}

point_transformation mrf_sampler::get_swap_transformation(size_t move_from_unit,
    size_t move_to_unit) {
  string from_label = unit_labels[move_from_unit];
  string to_label = unit_labels[move_to_unit];
  string swap_key = from_label + "-" + to_label;
  return swap_transformations[swap_key];
}

vector<pair<size_t, point_transformation_sequence> >
    mrf_sampler::get_all_adjusted_transformations(
        size_t unit_index, point_transformation base_transformation,
        bool swapping_enabled) {
  vector<pair<size_t, point_transformation_sequence> > all;
  if (swapping_enabled) {
    for (size_t alternate_index = 0; alternate_index < aligned_units.size();
         alternate_index++) {
      if (alternate_index != unit_index) {
        point_transformation_sequence adjusted_sequence;
        point_transformation swap_transformation =
            get_swap_transformation(alternate_index, unit_index);
        add_centroid_adjusted_transformation(alternate_index,
            base_transformation, &adjusted_sequence);
        adjusted_sequence.add(swap_transformation);
        all.push_back(std::make_pair(alternate_index, adjusted_sequence));
      }
    }
  } else {
    point_transformation_sequence adjusted_sequence;
    add_centroid_adjusted_transformation(unit_index,
          base_transformation, &adjusted_sequence);
    all.push_back(std::make_pair(unit_index, adjusted_sequence));
  }
  return all;
}

void mrf_sampler::add_centroid_adjusted_transformation(size_t unit_index,
    point_transformation transformation,
    point_transformation_sequence* sequence) {
  double centroid_x, centroid_y, centroid_z;
  aligned_units[unit_index].get_centroid(&centroid_x, &centroid_y, &centroid_z);
  sequence->add_translation(-centroid_x, -centroid_y, -centroid_z);
  sequence->add(transformation);
  sequence->add_translation(centroid_x, centroid_y, centroid_z);
}
