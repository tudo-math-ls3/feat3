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
* MPI_Dist_graph_create(...), (see MPI-2.2 standard, example 7.3 on page 256).
*
* COMMENT_HILMAR: This is only a very rough first version, which will be surely adapted to our needs...
*
* COMMENT_HILMAR: This graph structure is the most general one. When it comes to MPI communication, we surely have to
*   distinguish edge neighbours and diagonal neighbours.
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
  int const _num_nodes;

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
  * The neighbours of node i are stored in the subarray _edges[#_index[i]], ..., _edges[#_index[i+1]-1]. The order
  * of the neighbours within the subarrays is arbitrary.
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
  * Dimension: [total number of neighbours] = [number of edges] = [#_index[#_num_nodes]]
  */
  int* _edges;


public:

  /* *****************
  * member variables *
  *******************/

  /* *************************
  * constructor & destructor *
  ***************************/
  /// constructor
  Graph(
    int const num_nodes,
    int* index,
    int* edges
    )
    : _num_nodes(num_nodes),
      _index(index),
      _edges(edges)
  {
  }
  /// destructor
  ~Graph()
  {
    delete [] _index;
    _index = nullptr;
    delete [] _edges;
    _edges = nullptr;
  }

  /* ******************
  * getters & setters *
  ********************/
  /**
  * \brief getter for the number of nodes
  *
  * \return number of nodes #_num_nodes
  */
  inline int num_nodes() const
  {
    return _num_nodes;
  }

  /**
  * \brief getter for the index array
  *
  * \return pointer to the index array #_index
  */
  inline int* index() const
  {
    return _index;
  }

  /**
  * \brief getter for the edge array
  *
  * \return pointer to the edge array #_edges
  */
  inline int* edges() const
  {
    return _edges;
  }

  /* *****************
  * member functions *
  *******************/
  /// print the graph
  void print()
  {
    std::cout << "number of nodes: " << _num_nodes << std::endl;
    if (_num_nodes > 0)
    {
      std::cout << "node | degree | neighbours: " << std::endl;
      for(int i(0) ; i < _num_nodes ; ++i)
      {
        std::cout << i << " | " << _index[i+1] - _index[i];
        if (_index[i+1] - _index[i] > 0)
        {
          std::cout << " | " << _edges[_index[i]];
          for(int j(_index[i]+1) ; j < _index[i+1] ; ++j)
          {
            std::cout << ", " << _edges[j];
          }
        }
        std::cout << std::endl;
      }
    }
  }
}; // class Graph



/**
* \brief class providing a distributed graph data structure for defining connectivity of subdomains / matrix patches /
*        processes
*
* This data structure represents the part of the global graph this process is associated with. It basically stores the
* number of neighbours and the ranks of the neighbours. The data structure can be constructed from a global graph with
* the help of the MPI functions
*   MPI_Dist_graph_neighbors_count(...)
* to get the number of neighbours and
*   MPI_Dist_graph_neighbors(...)
* to get the ranks of the neighbours. With the help of this data structure the global MPI topology graph can be
* created via the function MPI_Dist_graph_create(...).
*
* COMMENT_HILMAR: This is only a very rough first version, which will be surely adapted to our needs...
*
* \author Hilmar Wobker
*/
class GraphDistributed
{

private:

  /* *****************
  * member variables *
  *******************/
  /// number of neighbours
  int _num_neighbours;

  /**
  * \brief ranks of the neighbours
  *
  * Dimension: [#_num_neighbours]
  */
  int* _neighbours;


public:

  /* *************************
  * constructor & destructor *
  ***************************/
  /// constructor
  GraphDistributed(
    int const num_neighbours,
    int* neighbours
    )
    : _num_neighbours(num_neighbours),
      _neighbours(neighbours)
  {
  }

  /// destructor
  ~GraphDistributed()
  {
    delete [] _neighbours;
    _neighbours = nullptr;
  }

  /* ******************
  * getters & setters *
  ********************/
  /**
  * \brief getter for the number of neighbours
  *
  * \return number of neighbours #_num_neighbours
  */
  inline int num_neighbours() const
  {
    return _num_neighbours;
  }

  /**
  * \brief getter for the neighbour array
  *
  * \return pointer to the neibhbour array #_neighbours
  */
  inline int* neighbours() const
  {
    return _neighbours;
  }

  /* *****************
  * member functions *
  *******************/
}; // class GraphDistributed

#endif // guard KERNEL_GRAPH_HPP
