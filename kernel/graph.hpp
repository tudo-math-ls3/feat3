#pragma once
#ifndef KERNEL_GRAPH_HPP
#define KERNEL_GRAPH_HPP 1

// includes, system
#include <iostream>
#include <stdlib.h>

// includes, Feast

/**
* \brief class providing a graph data structure for defining connectivity of subdomains / matrix patches / processes
*
* This data structure is very similar to that described in the MPI-2.2 standard (p. 250ff) for creating process
* topologies (for the difference, see description of member #_index).
*
* In the case this data structure is used for storing the MPI communication graph, this is a \em global representation,
* i.e. a process storing this graph knows the \em complete communication graph. (For a distributed representation, see
* data structure GraphDistributed.) To construct a corresponding MPI communicator, one can use the function
*   int MPI_Graph_create(MPI_Comm comm_old, int num_nodes, int *index, int *edges, int reorder, MPI_Comm *comm_graph)
* (see MPI-2.2 standard, p. 250). When passing the array #_index to this routine, remember to omit the first entry
* (see description of the array #_index).
*
* The data structure can also be used to create distributed graph data structures via the MPI routine
* MPI_DIST_GRAPH_CREATE(...), (see MPI-2.2 standard, example 7.3 on page 256).
*
* COMMENT_HILMAR: Maybe a clever procedure: coordinator process creates global graph structure, distributes it to
* all the other processes of the process group via MPI_Dist_graph_create(...). Then, each process creates its local
* graph structure by inquiring corresponding information. (MPI_Topo_test, MPI_Graph_dims_get, MPI_Graph_get,
* MPI_Graph_neighbors_count, MPI_Graph_neighbors, MPI_Dist_graph_neighbors_count, MPI_Dist_graph_neighbors,
*
* COMMENT_HILMAR: This is only a very rough first version, which will be surely adapted to our needs...
*
* \author Hilmar Wobker
*/
class Graph
{

private:

  /* *****************
  * member variables *
  *******************/
  /// number of nodes in the graph, which are numbered from 0 to num_nodes-1
  int _num_nodes;

  /**
  * \brief access information for the array edges
  *
  * The i-th entry of this array stores the total number of neighbours of the graph nodes 0, ..., i-1, where
  * index[0] = 0. So, the i-th entry gives the position in the array #_edges, in which the subarray storing the node
  * neighbours of node i starts. The degree (=number of neighbours) of node i is given by _index[i+1] - _index[i].
  * (Difference to the data structure described in the MPI-2.2 standard: There, the 0-th entry is omitted. However,
  * adding this entry eases array access since the first node does not have to be handled in a special way.)
  * For an example see #_edges.
  *
  * Dimension: [#_num_nodes+1]
  */
  int* _index;

  /**
  * \brief edges of the graph, represented as a list of node numbers
  *
  * The neighbours of node i are stored in the subarray _edges[#_index[i]], ..., _edges[#_index[i+1]-1].
  *
  * Example: domain consisting of 7 subdomains, each subdomain is a node in the graph. Edges represent neighbouring
  * subdomains, including diagonal neighbours.
  *   -----------------
  *   | 0 | 1 | 2 | 3 |
  *   -----------------
  *   | 4 | 5 |
  *   ---------
  *   | 6 |
  *   -----
  *   nodes     neighbours
  *   0         4,5,1
  *   1         0,4,5,2
  *   2         1,5,3
  *   3         2
  *   4         0,1,5,6
  *   5         0,1,2,4,6
  *   6         4,5
  * _num_nodes = 7
  * _index = [0,        3,           7,      10, 11,          15,          20,    22]
  * _edges = [4, 5, 1,  0, 4, 5, 2,  1, 5, 3, 2,  0, 1, 5, 6,  0, 1, 2, 4,  6, 4,  5]
  *
  * Dimension: [total number of neighbours] = [number of edges] = [#_index[i+1]]
  */
  int* _edges;


public:

  /* *****************
  * member variables *
  *******************/

  /* *************************
  * constructor & destructor *
  ***************************/

  /* *****************
  * member functions *
  *******************/
}; // class Graph

#endif // guard KERNEL_GRAPH_HPP
