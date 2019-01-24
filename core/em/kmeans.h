#ifndef _KMEANS_H_
#define _KMEANS_H_

#include <vector>

#include "../proteinstructure/MersenneTwister.h"

// T must implement "double* T.coordinates()" to extract the
// real-vector d-dimensional part of T that we care about for k-means.
// The caller must release the memory from the double* returned.
template<class T>
class kmeans_node
{
 private:
  // Dimensions corresponds to the size of coordinates.
  size_t dimensions;
  // Index from 0 to (N-1), from the N clusters managed by a
  // kmeans instance. This is used to determine if the clustering converged.
  size_t assigned_cluster;
  // This is a copy of the underlying real-vector contained in T.
  // Since T can have a non-array representation of them, to keep
  // the kmeans code generic we maintain a copy.
  double* coordinates;
  // A reference the the original instance. Useful if we just want to return
  // the objects that represent cluster centroids.
  T* base;
 public:
  kmeans_node();
  kmeans_node(const kmeans_node<T>& other);
  kmeans_node(T* in_base, size_t in_dimensions, double* in_coordinates);
  ~kmeans_node();
  size_t get_assigned_cluster();
  void set_assigned_cluster(size_t new_assignment);
  double* get_coordinates() const;
  T* get_base_object();
  // Euclidean distance.
  double distance(double* other_coordinates) const;
};

template<class T>
class kmeans_cluster
{
 private:
  // This is a pointer to all the kmeans_node instances stored in the parent
  // kmeans object. The actual elements that belong to the current
  // kmeans_cluster are stored as indices in nodes_in_cluster
  vector<kmeans_node<T> >* node_pointers;
  vector<size_t> nodes_in_cluster;
  // This is the average of the coordinates. It is not necessarily the exact
  // coordinate of the a node.
  double* centroid_coordinates;
  size_t dimensions;
  // Index within node_pointers to the node that is the closest
  // to centroid_coordinates
  size_t centroid_node_index;
  
 public:
  kmeans_cluster();
  kmeans_cluster(const kmeans_cluster<T>& other);
  // Maintains a reference to the main node collection and allocates
  // the memory needed for the cluster centroid, which must be provided by
  // selecting one randomly (performed by a kmeans instance).
  kmeans_cluster(vector<kmeans_node<T> >* in_node_pointers,
                 size_t in_dimensions, size_t random_centroid_index);
  // Just release the centroid memory.
  ~kmeans_cluster();
  // Returns the average coordinates, not the actual node.
  double* get_centroid_coordinates();
  // This returns the index in the original collection that represents
  // the node that is closest to the centroid coordinates.
  size_t get_centroid_index();
  size_t get_size();
  // This is based on the "nodes" average and overwrites both
  // centroid_coordinates and centroid_node_index
  void recalculate_centroid();
  // Should be called before the k-means assignment step. Every node will be
  // reassigned based on the current values of the centroids.
  void clear_nodes_assigned();
  // Adds to nodes_in_cluster.
  void add_node(size_t node_index);
  // Returns the average distance between the node represented by node_index
  // and all nodes in this cluster. If node index is already included in this
  // cluster it will be ignored. This value is used to create silhouette
  // coefficients.
  double calculate_average_node_dissimilarity(size_t node_index);
  // Max distance between the centroid and a single node in the cluster
  double calculate_radius();
  // Max distance between any two nodes in the cluster
  double calculate_diameter();
};

template<class T>
class kmeans
{
 private:
  // All nodes clustered are expected to have these many dimensions.
  size_t dimensions;
  vector<kmeans_node<T> > nodes;
  vector<kmeans_cluster<T> > clusters;
  // To randomly initialize the cluster centers
  MTRand random_generator;
  // Auxiliary method use by calculate_silhouette_coefficient that
  // computes the silhouette for just one node.
  double calculate_silhouette_for_node(size_t node_index);
  
 public:
  kmeans();
  kmeans(size_t in_dimensions);
  // Simply add the pointer to the nodes.
  void add_node(T* new_node);
  // Initial k centroids and k empty clusters, randomly from the current nodes.
  void initialize(size_t k);
  // TODO: set iteration limit for convergence purposes?
  void cluster();
  // Should only be called after cluster()
  vector<T*> get_centroids();
  vector<size_t> get_cluster_sizes();
  // Calculates the overall silhouette coefficient for all data points
  // clustered.
  // http://en.wikipedia.org/wiki/Silhouette_(clustering)
  double calculate_silhouette_coefficient();
  double calculate_average_cluster_radius();
  double calculate_average_cluster_diameter();
};

#endif
