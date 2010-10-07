#pragma once
#ifndef KERNEL_BASE_MESH_HPP
#define KERNEL_BASE_MESH_HPP 1

// includes, system
#include <iostream>
#include <stdlib.h>

// includes, Feast
#include <kernel/graph.hpp>

/**
* \brief class describing a base mesh
*
* \author Hilmar Wobker
*/
class BaseMesh
{

private:

  /* *****************
  * member variables *
  *******************/
  /// graph describing the connectivity of the base mesh
  Graph* _graph;


public:

  /* *****************
  * member variables *
  *******************/

  /* *************************
  * constructor & destructor *
  ***************************/
  BaseMesh()
    : _graph(nullptr)
  {
  }

  ~BaseMesh()
  {
    if (_graph != nullptr)
    {
      delete _graph;
      _graph = nullptr;
    }
  }

  /* ******************
  * getters & setters *
  ********************/
  /**
  * \brief getter for the graph
  *
  * \return Graph pointer #_graph
  */
  inline Graph* graph() const
  {
    return _graph;
  }

  /* *****************
  * member functions *
  *******************/
  void read_mesh()
  {
    // read mesh
    // ...

    // not implemented yet, so manually create a graph describing some base mesh
    // subdomain  num_neigh   neigh             num_diag    diag
    //    0          4        3,4,5,1              2        4,5
    //    1          6        0,3,4,5,6,2          3        3,4,6
    //    2          4        1,5,6,7              2        5,7
    //    3          6        8,9,4,5,1,0          3        9,5,1
    //    4          6        5,1,0,3,8,9          3        1,0,8
    //    5          6        6,2,1,0,3,4          3        2,0,3
    //    6          4        7,2,1,5              1        1
    //    7          3        12,2,6               1        2
    //    8          6        13,14,10,9,4,3       3        14,10,4
    //    9          6        4,3,8,13,14,10       3        3,13,14
    //   10          6        9,8,13,14,15,11      3        8,13,15
    //   11          4        10,14,15,12          1        14
    //   12          3        11,15,7              1        15
    //   13          4        14,10,9,8            2        10,9
    //   14          6        15,11,10,9,8,13      3        11,9,8
    //   15          4        12,11,10,14          2        12,10
    int const num_base_cells = 16;
    int* index = new int[num_base_cells+1];
    index[0]  =  0;
    index[1]  =  4;
    index[2]  = 10;
    index[3]  = 14;
    index[4]  = 20;
    index[5]  = 26;
    index[6]  = 32;
    index[7]  = 36;
    index[8]  = 39;
    index[9]  = 45;
    index[10] = 51;
    index[11] = 57;
    index[12] = 61;
    index[13] = 64;
    index[14] = 68;
    index[15] = 74;
    index[16] = 78;

    int* edges = new int[index[num_base_cells]];
    // neighbours of subdomain 0
    edges[0]  =  3;  edges[1] =  4;  edges[2] =  5;  edges[3] =  1;
    // neighbours of subdomain 1
    edges[4]  =  0;  edges[5] =  3;  edges[6] =  4;  edges[7] =  5;  edges[8] =  6;  edges[9] =  2;
    // neighbours of subdomain 2
    edges[10] =  1; edges[11] =  5; edges[12] =  6; edges[13] =  7;
    // neighbours of subdomain 3
    edges[14] =  8; edges[15] =  9; edges[16] =  4; edges[17] =  5;  edges[18] =  1; edges[19] =  0;
    // neighbours of subdomain 4
    edges[20] =  5; edges[21] =  1; edges[22] =  0; edges[23] =  3;  edges[24] =  8; edges[25] =  9;
    // neighbours of subdomain 5
    edges[26] =  6; edges[27] =  2; edges[28] =  1; edges[29] =  0;  edges[30] =  3; edges[31] =  4;
    // neighbours of subdomain 6
    edges[32] =  7; edges[33] =  2; edges[34] =  1; edges[35] =  5;
    // neighbours of subdomain 7
    edges[36] = 12; edges[37] =  2; edges[38] =  6;
    // neighbours of subdomain 8
    edges[39] = 13; edges[40] = 14; edges[41] = 10; edges[42] =  9;  edges[43] =  4; edges[44] =  3;
    // neighbours of subdomain 9
    edges[45] =  4; edges[46] =  3; edges[47] =  8; edges[48] = 13;  edges[49] = 14; edges[50] = 10;
    // neighbours of subdomain 10
    edges[51] =  9; edges[52] =  8; edges[53] = 13; edges[54] = 14;  edges[55] = 15; edges[56] = 11;
    // neighbours of subdomain 11
    edges[57] = 10; edges[58] = 14; edges[59] = 15; edges[60] = 12;
    // neighbours of subdomain 12
    edges[61] = 11; edges[62] = 15; edges[63] =  7;
    // neighbours of subdomain 13
    edges[64] = 14; edges[65] = 10; edges[66] =  9; edges[67] =  8;
    // neighbours of subdomain 14
    edges[68] = 15; edges[69] = 11; edges[70] = 10; edges[71] =  9;  edges[72] =  8; edges[73] = 13;
    // neighbours of subdomain 15
    edges[74] = 12; edges[75] = 11; edges[76] = 10; edges[77] = 14;

    if (_graph != nullptr)
    {
      delete _graph;
      _graph = nullptr;
    }
    // create graph object
    // temporarily, do not distinguish edge neighbours and diagonal neighbours
    _graph = new Graph(num_base_cells, index, edges);
    _graph->print();
  }
}; // class BaseMesh

#endif // guard KERNEL_BASE_MESH_HPP
