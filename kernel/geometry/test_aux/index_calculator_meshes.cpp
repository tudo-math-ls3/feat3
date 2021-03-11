// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/test_aux/index_calculator_meshes.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>

namespace FEAT
{
  namespace Geometry
  {
    namespace TestAux
    {
      void validate_hypercube_edge_index_set(const HEdgeIndexTree& tree)
      {
        // vertices-at-edge array
        static const Index v_e[] =
        {
          0, 1, // v0
          1, 3,
          10, 8,
          42, // v1
          2, 4,
          11, 9,
          42, //v2
          3, 3,
          5, 5,
          12, 10,
          42, //v3
          4, 4,
          6, 6,
          13, 11,
          42, //v4
          7, 7,
          14, 12,
          42, //v5
          8, 6,
          15, 13,
          42, //v6
          9, 7,
          16, 14,
          42, //v7
          17, 15,
          42, //v8
          18, 9,
          19, 11,
          28, 16,
          42, //v9
          20, 12,
          29, 17,
          42, //v10
          21, 11,
          23, 13,
          42, //v11
          22, 12,
          24, 14,
          30, 18,
          42, //v12
          25, 15,
          31, 19,
          42, //v13
          26, 14,
          42, //v14
          27, 15,
          42, //v15
          42, //v16
          32, 17,
          33, 18,
          42, //v17
          34, 19,
          42, //v18
          35, 19,
          42, //v19
          42
        };

        // checking
        Index count = 0;
        Index current_value;
        Index vector_size = 19;
        int num_indices = tree.get_num_indices();

        for(Index i(0); i < vector_size; ++i)
        {
          Index current_set_size = tree.get_set_size(i);
          for(Index j(0); j < current_set_size; ++j)
          {
            for(int k(0); k < num_indices; ++k)
            {
              current_value = tree.get_index(i, j, k);

              // if there is a wrong value
              if(current_value != v_e[count])
                throw String("IndexTree hypercube-edge-mesh parsing failure");
              ++count;
            }
          }

          // if something is missing
          if (v_e[count] !=42)
            throw String("IndexTree hypercube-edge-mesh parsing failure");
          count++;
        }

      } // validate_hypercube_edge_index_set

      void validate_hypercube_quad_index_set(const QuadIndexTree& tree)
      {

        // vertices-at-quad array
        static const Index v_q[] =
        {
          0, 1, 3, 4,
          3, 1, 8, 9,
          4, 3, 8, 11,
          42,
          5, 4, 9, 12,
          42,
          1, 3, 5, 6,
          6, 3, 10, 11,
          8, 5, 10, 13,
          42,
          2, 4, 6, 7,
          7, 4, 11, 12,
          9, 6, 11, 14,
          42,
          10, 7, 12, 15,
          42,
          11, 6, 13, 14,
          42,
          12, 7, 14, 15,
          42,
          42,
          13, 9, 11, 12,
          16, 9, 16, 17,
          17, 11, 16, 18,
          42,
          18, 12, 17, 19,
          42,
          14, 11, 13, 14,
          42,
          15, 12, 14, 15,
          19, 12, 18, 19,
          42,
          42,
          42,
          42,
          42,
          20, 17, 18, 19,
          42,
          42,
          42
        };

        // checking
        Index count = 0;
        Index current_value;
        Index vector_size = 19;
        int num_indices = tree.get_num_indices();

        for(Index i(0); i < vector_size; ++i)
        {
          Index current_set_size = tree.get_set_size(i);
          for(Index j(0); j < current_set_size; ++j)
          {
            for(int k(0); k < num_indices; ++k)
            {
              current_value = tree.get_index(i, j, k);

              // if there is a wrong value
              if(current_value != v_q[count])
                throw String("IndexTree hypercube-quad-mesh parsing failure");
              ++count;
            }
          }

          // if something is missing
          if (v_q[count] !=42)
            throw String("IndexTree hypercube-quad-mesh parsing failure");
          count++;
        }

      } // validate_hypercube_quad_index_set

      void validate_simplex_edge_index_set(const SEdgeIndexTree& tree)
      {

        // vertices-at-edge array
        static const Index v_e[] =
        {
          7, 1,
          9, 2,
          4, 3,
          1, 4,
          42,
          8, 2,
          5, 3,
          2, 4,
          42,
          6, 3,
          3, 4,
          42,
          0, 4,
          42,
          42
        };

        //checking
        Index count = 0;
        Index current_value;
        Index vector_size = 4;
        int num_indices = tree.get_num_indices();

        for(Index i(0); i < vector_size; ++i)
        {
          Index current_set_size = tree.get_set_size(i);
          for(Index j(0); j < current_set_size; ++j)
          {
            for(int k(0); k < num_indices; ++k)
            {
              current_value = tree.get_index(i, j, k);

              // if there is a wrong value
              if(current_value != v_e[count])
                throw String("IndexTreesimplex-edge-mesh parsing failure");
              ++count;
            }
          }

          // if something is missing
          if (v_e[count] !=42)
            throw String("IndexTree simplex-edge-mesh parsing failure");
          count++;
        }

      } // validate_simplex_edge_index_set

      void validate_simplex_triangle_index_set(const TriaIndexTree& tree)
      {

        // vertices-at-triangle array
        static const Index v_t[] =
        {
          2, 1, 3,
          6, 1, 4,
          3, 2, 3,
          5, 2, 4,
          8, 3, 4,
          42,
          1, 2, 3,
          4, 2, 4,
          0, 3, 4,
          42,
          7, 3, 4,
          42,
          42
        };

        //checking
        Index count = 0;
        Index current_value;
        Index vector_size = 4;
        int num_indices = tree.get_num_indices();

        for(Index i(0); i < vector_size; ++i)
        {
          Index current_set_size = tree.get_set_size(i);
          for(Index j(0); j < current_set_size; ++j)
          {
            for(int k(0); k < num_indices; ++k)
            {
              current_value = tree.get_index(i, j, k);

              // if there is a wrong value
              if(current_value != v_t[count])
                throw String("IndexTree simplex-triangle-mesh parsing failure");
              ++count;
            }
          }

          // if there is something missing
          if (v_t[count] !=42)
            throw String("IndexTree simplex-triangle-mesh parsing failure");
          count++;
        }

      } // validate_simplex_triangle_index_set

    } // namespace TestAux
  } // namespace Geometry
} // namespace FEAT
