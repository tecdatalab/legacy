#include "kmeans.h"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <limits>


// kmeans_node methods
template<class T>
kmeans_node<T>::kmeans_node() : dimensions(0), assigned_cluster(0),
    coordinates(NULL), base(NULL) {
}

template<class T>
kmeans_node<T>::kmeans_node(T* in_base, size_t in_dimensions,
    double* in_coordinates) : dimensions(in_dimensions), base(in_base){
  coordinates = new double[in_dimensions];
  memcpy(coordinates, in_coordinates, sizeof(double) * in_dimensions);
  // Set the assigned cluster to zero by default since we at least expect k=1.
  assigned_cluster = 0;
}

template<class T>
kmeans_node<T>::kmeans_node(const kmeans_node<T>& other) :
    dimensions(other.dimensions), assigned_cluster(other.assigned_cluster),
    base(other.base) {
  coordinates = new double[other.dimensions];
  memcpy(coordinates, other.coordinates, sizeof(double) * other.dimensions);
}

template<class T>
kmeans_node<T>::~kmeans_node() {
  if(coordinates != NULL) {
    delete [] coordinates;
  }
}

template<class T>
size_t kmeans_node<T>::get_assigned_cluster() {
  return assigned_cluster;
}

template<class T>
void kmeans_node<T>::set_assigned_cluster(size_t new_assignment) {
  assigned_cluster = new_assignment;
}

template<class T>
double* kmeans_node<T>::get_coordinates() const {
  return coordinates;
}

template<class T>
T* kmeans_node<T>::get_base_object() {
  return base;
}

template<class T>
double kmeans_node<T>::distance(double* other_coordinates) const {
  double sum = 0;
  for(size_t coordinate_index = 0; coordinate_index < this->dimensions;
      coordinate_index++) {
    double difference = this->coordinates[coordinate_index]
                        - other_coordinates[coordinate_index];
    sum += difference * difference;
  }
  return sqrt(sum);
}

// kmeans_cluster methods
template<class T>
kmeans_cluster<T>::kmeans_cluster() : node_pointers(NULL),
    centroid_coordinates(NULL), dimensions(0), centroid_node_index(0) {
}

template<class T>
kmeans_cluster<T>::kmeans_cluster(const kmeans_cluster<T>& other) :
    node_pointers(other.node_pointers),
    nodes_in_cluster(other.nodes_in_cluster), dimensions(other.dimensions),
    centroid_node_index(other.centroid_node_index) {
  centroid_coordinates = new double[other.dimensions];
  memcpy(centroid_coordinates, other.centroid_coordinates,
         sizeof(double) * other.dimensions);
}

template<class T>
kmeans_cluster<T>::kmeans_cluster(vector<kmeans_node<T> >* in_node_pointers,
    size_t in_dimensions, size_t random_centroid_index) :
        node_pointers(in_node_pointers), dimensions(in_dimensions),
        centroid_node_index(random_centroid_index) {
  centroid_coordinates = new double[in_dimensions];
  memcpy(centroid_coordinates,
         in_node_pointers->at(random_centroid_index).get_coordinates(),
         sizeof(double) * in_dimensions);
}

template<class T>
kmeans_cluster<T>::~kmeans_cluster() {
  if(centroid_coordinates != NULL) {
    delete [] centroid_coordinates;
  }
}

template<class T>
double* kmeans_cluster<T>::get_centroid_coordinates() {
  return centroid_coordinates;
}

template<class T>
size_t kmeans_cluster<T>::get_centroid_index() {
  return centroid_node_index;
}

template<class T>
size_t kmeans_cluster<T>::get_size() {
  return nodes_in_cluster.size();
}

template<class T>
void kmeans_cluster<T>::recalculate_centroid() {
  memset(centroid_coordinates, 0, sizeof(double) * dimensions);
  for (size_t in_cluster_index = 0; in_cluster_index < nodes_in_cluster.size();
       in_cluster_index++) {
    size_t overall_node_index = nodes_in_cluster[in_cluster_index];
    double* node_coordinates =
      node_pointers->at(overall_node_index).get_coordinates();
    // Add the coordinates to the cumulative we have so far.
    for (size_t dimension = 0; dimension < this->dimensions; dimension++) {
      centroid_coordinates[dimension] += node_coordinates[dimension];
    }  
  }
  // Divide by the total number of elements in the cluster to have the final
  // average value.
  for (size_t dimension = 0; dimension < this->dimensions; dimension++) {
    centroid_coordinates[dimension] /= nodes_in_cluster.size();
  }

  // Finally, with the recomputed centroid, find the node in the cluster that
  // is closest to it and assign it as the node centroid.
  double lowest_distance = std::numeric_limits<double>::max(); 
  for (size_t in_cluster_index = 0; in_cluster_index < nodes_in_cluster.size();
       in_cluster_index++) {
    size_t overall_node_index = nodes_in_cluster[in_cluster_index];
    double current_distance =
      node_pointers->at(overall_node_index).distance(centroid_coordinates);
    if(current_distance < lowest_distance) {
      centroid_node_index = overall_node_index;
      lowest_distance = current_distance;
    }
  } 
}

template<class T>
void kmeans_cluster<T>::clear_nodes_assigned() {
  nodes_in_cluster.clear();
}

template<class T>
void kmeans_cluster<T>::add_node(size_t node_index) {
  nodes_in_cluster.push_back(node_index);
}

template<class T>
double kmeans_cluster<T>::calculate_average_node_dissimilarity(
    size_t node_index) {
  // If it's found, the sum is divided by size-1 instead of size
  bool is_node_in_this_cluster = false;
  double sum = 0.0f;
  double* node_coordinates = node_pointers->at(node_index).get_coordinates();
  for(size_t current_index = 0; current_index < nodes_in_cluster.size();
      current_index++) {
    size_t current_node = nodes_in_cluster[current_index];
    if(current_node == node_index) {
      is_node_in_this_cluster = true;
    } else {
      sum += node_pointers->at(current_node).distance(node_coordinates);
    }
  }
  if(is_node_in_this_cluster) {
    return sum / static_cast<double>(nodes_in_cluster.size() - 1);
  } else {
    return sum / static_cast<double>(nodes_in_cluster.size());
  }
}

template<class T>
double kmeans_cluster<T>::calculate_radius() {
  double radius = 0.0f;
  for(size_t node_index = 0; node_index < nodes_in_cluster.size();
      node_index++) {
    const kmeans_node<T> &current_node =
        node_pointers->at(nodes_in_cluster[node_index]);
    double distance = current_node.distance(centroid_coordinates);
    if(distance > radius) {
      radius = distance;
    }
  }
  return radius; 
}

template<class T>
double kmeans_cluster<T>::calculate_diameter() {
  double diameter = 0.0f;
  for(size_t node1_index = 0; node1_index < nodes_in_cluster.size();
      node1_index++) {
    for(size_t node2_index = node1_index + 1;
        node2_index < nodes_in_cluster.size(); node2_index++) {
      const kmeans_node<T> &node1 =
          node_pointers->at(nodes_in_cluster[node1_index]);
      const kmeans_node<T> &node2 =
          node_pointers->at(nodes_in_cluster[node2_index]);
      double distance = node1.distance(node2.get_coordinates());
      if(distance > diameter) {
        diameter = distance;
      }
    }
  }
  return diameter;
}

// kmeans methods

template<class T>
double kmeans<T>::calculate_silhouette_for_node(size_t node_index) {
  size_t node_cluster = nodes[node_index].get_assigned_cluster();
  // a and b come from the terminology used in the wikipedia page
  // http://en.wikipedia.org/wiki/Silhouette_(clustering)
  double a =
      clusters[node_cluster].calculate_average_node_dissimilarity(node_index);
  double b = std::numeric_limits<double>::max();
  for(size_t cluster_index = 0; cluster_index < clusters.size();
      cluster_index++) {
    if(cluster_index == node_cluster) {
      continue;
      // because this is "double a" assigned before
    } else {
      double possible_b =
          clusters[cluster_index].calculate_average_node_dissimilarity(
              node_index);
      if(possible_b < b) {
        b = possible_b;
      }
    }
  }
  // return s = (b - a) / max(a,b)
  if(a > b) {
    return (b - a) / a;
  } else {
    return (b - a) / b;
  }
}

template<class T>
kmeans<T>::kmeans() : dimensions(0) {
}

template<class T>
kmeans<T>::kmeans(size_t in_dimensions) : dimensions(in_dimensions) {
}

template<class T>
void kmeans<T>::add_node(T* new_node) {
  // The memory needs to be released. The kmeans_node constructor
  // allocates internal memory in the instance and makes a copy, thus
  // the memory can be freed at the end of this method.
  double* node_coordinates = new_node->coordinates();
  kmeans_node<T> node(new_node, this->dimensions, node_coordinates);
  nodes.push_back(node); 
  delete [] node_coordinates;
}

template<class T>
void kmeans<T>::initialize(size_t k) {
  // First generate k random centroids.
  vector<size_t> random_centroids;
  while(random_centroids.size() < k) {
    size_t centroid_index =
        (size_t)floor(random_generator.randExc(nodes.size()));
    // Don't add an already generated centroid index. Just linear search
    // is fine most of the time, so we keep the implementation simple.
    if(std::find(random_centroids.begin(),
                 random_centroids.end(),
                 centroid_index) == random_centroids.end()) {
      random_centroids.push_back(centroid_index);
    }
  }
  // Empty the clusters and create new ones with the random centroids.
  clusters.clear();
  for(size_t cluster_counter = 0; cluster_counter < k; cluster_counter++) {
    kmeans_cluster<T> new_cluster(&nodes, this->dimensions,
                                  random_centroids[cluster_counter]);
    clusters.push_back(new_cluster); 
  } 
}

template<class T>
void kmeans<T>::cluster() {
  // Repeat the process while the assignments change.
  bool assignments_changed;
  do {
    assignments_changed = false;
    // 1. Clear all previous node assignments.
    for(size_t cluster_index = 0; cluster_index < clusters.size();
        cluster_index++) {
      clusters[cluster_index].clear_nodes_assigned();
    }
    // 2. Find shortest distance to centroid for each node and assign.
    for(size_t node_index = 0; node_index < nodes.size(); node_index++) {
      double lowest_distance = std::numeric_limits<double>::max();
      // Default cluster assignment should change immediately because everything
      // will be lower than double_max
      size_t current_cluster_assignment = 0;
      for(size_t cluster_index = 0; cluster_index < clusters.size();
          cluster_index++) {
        double current_distance =
            nodes[node_index].distance(
                clusters[cluster_index].get_centroid_coordinates());
        if(current_distance < lowest_distance) {
          current_cluster_assignment = cluster_index;
          lowest_distance = current_distance;
        }
      }
      clusters[current_cluster_assignment].add_node(node_index);
      // Check if the node assignment changed from the previous iteration.
      // Note that if the flag was previously true, we want to keep it "on",
      // thus the ||=.
      assignments_changed = assignments_changed ||
                            (nodes[node_index].get_assigned_cluster() !=
                             current_cluster_assignment);
      nodes[node_index].set_assigned_cluster(current_cluster_assignment);
    }
    // 3. Recompute the centroids.
    for(size_t cluster_index = 0; cluster_index < clusters.size();
        cluster_index++) {
      clusters[cluster_index].recalculate_centroid();
    }
  } while(assignments_changed);
}

template<class T>
vector<T*> kmeans<T>::get_centroids() {
  vector<T*> centroids;
  for(size_t cluster_index = 0; cluster_index < clusters.size();
      cluster_index++) {
    size_t centroid_index = clusters[cluster_index].get_centroid_index();
    centroids.push_back(nodes[centroid_index].get_base_object());
  }
  return centroids;
}

template<class T>
vector<size_t> kmeans<T>::get_cluster_sizes() {
  vector<size_t> sizes;
  for(size_t cluster_index = 0; cluster_index < clusters.size();
      cluster_index++) {
    sizes.push_back(clusters[cluster_index].get_size());
  }
  return sizes;
}

template<class T>
double kmeans<T>::calculate_silhouette_coefficient() {
  double sum = 0.0f;
  for(size_t node_index = 0; node_index < nodes.size(); node_index++) {
    sum += calculate_silhouette_for_node(node_index);
  }
  return sum / static_cast<double>(nodes.size());
}
template<class T>
double kmeans<T>::calculate_average_cluster_radius() {
  double sum = 0.0f;
  for(size_t cluster_index = 0; cluster_index < clusters.size();
      cluster_index++) {
    sum += clusters[cluster_index].calculate_radius();
  }
  return sum / static_cast<double>(clusters.size());
}

template<class T>
double kmeans<T>::calculate_average_cluster_diameter() {
  double sum = 0.0f;
  for(size_t cluster_index = 0; cluster_index < clusters.size();
      cluster_index++) {
    sum += clusters[cluster_index].calculate_diameter();
  }
  return sum / static_cast<double>(clusters.size());
}
