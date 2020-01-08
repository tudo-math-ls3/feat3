// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_STANDARD_INDEX_REFINER_HPP
#define KERNEL_GEOMETRY_INTERN_STANDARD_INDEX_REFINER_HPP 1

// includes, FEAT
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/intern/entity_counter.hpp>
#include <kernel/geometry/intern/standard_refinement_traits.hpp>
#include <kernel/geometry/intern/sub_index_mapping.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Shape_,
        int cell_dim_,
        int face_dim_>
      struct StandardIndexRefiner;

      /**
       * \brief StandardIndexRefiner implementation for Simplex<1>: Vertices-at-Edges
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Simplex<1>, 1, 0>
      {
        static constexpr int shape_dim = 1;
        static constexpr int cell_dim = 1; // second template parameter
        static constexpr int face_dim = 0;  // third template parameter

        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioe = index_offsets[1];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set type; an edge has 2 vertices, so we need an IndexSet<2>
          typedef IndexSet<2> IndexSetTypeEV;

          // typedef index vector reference type
          typedef IndexSetTypeEV::IndexTupleType IndexTupleTypeEV;

          // fetch the vertices-at-edge index set
          const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();

          // fetch number of coarse mesh edges
          const Index num_edges = index_set_e_v.get_num_entities();

          // loop over all coarse mesh edges
          for(Index i(0); i < num_edges; ++i)
          {
            // fetch coarse mesh vertices-at-edge index vector
            const IndexTupleTypeEV& e_v = index_set_e_v[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexTupleType& idx0 = index_set_out[offset + 2*i + 0];
            IndexTupleType& idx1 = index_set_out[offset + 2*i + 1];

            // calculate fine edge-vertex indices
            idx0[0] = iov + e_v[0];
            idx0[1] = ioe + i;
            idx1[0] = ioe + i;
            idx1[1] = iov + e_v[1];
          }

          // return fine edge count
          return 2*num_edges;
        }
      }; // StandardIndexRefiner<Simplex<1>,1,0>

      /* ************************************************************************************* */
      /* ************************************************************************************* */
      /* ************************************************************************************* */

      /**
       * \brief StandardIndexRefiner implementation for Simplex<2>: Vertices-at-Edges
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Simplex<2>, 1, 0>
      {
        static constexpr int shape_dim = 2;
        static constexpr int cell_dim = 1; // second template parameter
        static constexpr int face_dim = 0;  // third template parameter

        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index ioe = index_offsets[1];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set type; a triangle has 3 edges, so we need an IndexSet<3>
          typedef IndexSet<3> IndexSetTypeSE;

          // typedef index vector type
          typedef IndexSetTypeSE::IndexTupleType IndexTupleTypeSE;

          // fetch the edges-at-simplex index set
          const IndexSetTypeSE& index_set_s_e = index_set_holder_in.get_index_set<2,1>();

          // fetch number of coarse mesh simplices
          const Index num_simps = index_set_s_e.get_num_entities();

          /*
              v_2
              |\
              | \
              |  \
              |   \
              |    \
              e_1->-e_0
              |\    |\
              | \   | \
              |  ^  v  \
              |   \ |   \
              |    \|     \
              v_0---e_2----v_1
          */

          // loop over all coarse mesh simplices
          for(Index i(0); i < num_simps; ++i)
          {
            // fetch coarse mesh edges-at-triangle index vector
            const IndexTupleTypeSE& s_e = index_set_s_e[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexTupleType& e_0 = index_set_out[offset + 3*i + 0];
            IndexTupleType& e_1 = index_set_out[offset + 3*i + 1];
            IndexTupleType& e_2 = index_set_out[offset + 3*i + 2];

            e_0[0] = ioe + s_e[2]; // e_2
            e_0[1] = ioe + s_e[1]; // e_1

            e_1[0] = ioe + s_e[0]; // e_0
            e_1[1] = ioe + s_e[2]; // e_2

            e_2[0] = ioe + s_e[1]; // e_1
            e_2[1] = ioe + s_e[0]; // e_0

          }
          // return fine edge count
          return 3*num_simps;
        }
      }; // StandardIndexRefiner<Simplex<2>,1,0>

      /**
       * \brief StandardIndexRefiner implementation for Simplex<2>: Vertices-at-Tris
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Simplex<2>, 2, 0>
      {
        static constexpr int shape_dim = 2;
        static constexpr int cell_dim = 2; // second template parameter
        static constexpr int face_dim = 0;  // third template parameter

        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioe = index_offsets[1];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types; a triangle has 3 vertices and 3 edges, so both need an IndexSet<3>
          typedef IndexSet<3> IndexSetTypeSV;
          typedef IndexSet<3> IndexSetTypeSE;

          // typedef index vector type
          typedef IndexSetTypeSV::IndexTupleType IndexTupleTypeSV;
          typedef IndexSetTypeSE::IndexTupleType IndexTupleTypeSE;

          // fetch the simplex-vertex index set
          const IndexSetTypeSV& index_set_s_v = index_set_holder_in.get_index_set<2,0>();
          const IndexSetTypeSE& index_set_s_e = index_set_holder_in.get_index_set<2,1>();

          // fetch number of coarse mesh simplices
          const Index num_simps = index_set_s_v.get_num_entities();

          /*
              v_2
              |\
              | \
              |s \
              |2  \
              |    \
              e_1->-e_0
              |\  s |\
              | \ 3 | \
              |  ^  v  \
              |s  \ | s \
              |0   \| 1  \
              v_0---e_2----v_1
          */

          // loop over all coarse mesh simplices
          for(Index i(0); i < num_simps; ++i)
          {
            // fetch coarse mesh vertices-at-triangle and edges-at-triangle index vectors
            const IndexTupleTypeSV& s_v = index_set_s_v[i];
            const IndexTupleTypeSE& s_e = index_set_s_e[i];

            // fetch fine mesh vertices-at-triangle index vectors
            IndexTupleType& s_0 = index_set_out[offset + 4*i + 0];
            IndexTupleType& s_1 = index_set_out[offset + 4*i + 1];
            IndexTupleType& s_2 = index_set_out[offset + 4*i + 2];
            IndexTupleType& s_3 = index_set_out[offset + 4*i + 3];

            // calculate fine simplex-vertex indices
            s_0[0] = iov + s_v[0]; // v_0
            s_0[1] = ioe + s_e[2]; // e_2
            s_0[2] = ioe + s_e[1]; // e_1

            s_1[0] = ioe + s_e[2]; // e_2
            s_1[1] = iov + s_v[1]; // v_1
            s_1[2] = ioe + s_e[0]; // e_0

            s_2[0] = ioe + s_e[1]; // e_1
            s_2[1] = ioe + s_e[0]; // e_0
            s_2[2] = iov + s_v[2]; // s_2

            s_3[0] = ioe + s_e[0]; // e_0
            s_3[1] = ioe + s_e[1]; // e_1
            s_3[2] = ioe + s_e[2]; // e_2
          }
          // return fine simplex count
          return 4*num_simps;
        }
      }; // StandardIndexRefiner<Simplex<2>,2,0>

      /**
       * \brief StandardIndexRefiner implementation for Simplex<2>: Edges-at-Tris
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Simplex<2>, 2, 1>
      {
        static constexpr int shape_dim = 2;
        static constexpr int cell_dim = 2; // second template parameter
        static constexpr int face_dim = 1;  // third template parameter

        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch edge index offsets
          const Index ioe = index_offsets[1];
          const Index ios = index_offsets[2];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types;
          typedef IndexSet<2> IndexSetTypeEV; // an edge has 2 vertices
          typedef IndexSet<3> IndexSetTypeSV; // a simplex has 3 vertices
          typedef IndexSet<3> IndexSetTypeSE; // a simplex has 3 edges

          // typedef index vector type
          typedef IndexSetTypeSV::IndexTupleType IndexTupleTypeSV;
          typedef IndexSetTypeSE::IndexTupleType IndexTupleTypeSE;

          // fetch the index sets
          const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();
          const IndexSetTypeSV& index_set_s_v = index_set_holder_in.get_index_set<2,0>();
          const IndexSetTypeSE& index_set_s_e = index_set_holder_in.get_index_set<2,1>();

          // fetch number of coarse mesh simplices
          const Index num_simps = index_set_s_v.get_num_entities();

          // typedef the sub-index mapping
          typedef Intern::SubIndexMapping<ShapeType, face_dim, 0> SubIndexMappingType;

          /*
              v_2
              |\
              | \
              |s \
              |2  ^
              |    \
              e_1->-e_0
              |\  s |\
              v \ 3 | \
              |  ^  v  \
              |s  \ | s \
              |0   \| 1  \
              v_0---e_2-->-v_1
          */

          // loop over all coarse mesh simplices
          for(Index i(0); i < num_simps; ++i)
          {
            // fetch coarse mesh index vectors
            const IndexTupleTypeSV& s_v = index_set_s_v[i];
            const IndexTupleTypeSE& s_e = index_set_s_e[i];

            // create a sub-index mapping object
            SubIndexMappingType sim(s_v, s_e, index_set_e_v);

            // fetch fine mesh edges-at-simplex index vectors
            IndexTupleType& s_0 = index_set_out[offset + 4*i + 0];
            IndexTupleType& s_1 = index_set_out[offset + 4*i + 1];
            IndexTupleType& s_2 = index_set_out[offset + 4*i + 2];
            IndexTupleType& s_3 = index_set_out[offset + 4*i + 3];

            // calculate fine edges-at-simplex indices
            s_0[0] = ios + 3*i + 0;
            s_0[1] = ioe + 2*s_e[1] + sim.map(1, 1);  // e_1_1
            s_0[2] = ioe + 2*s_e[2] + sim.map(2, 0);  // e_2_0

            s_1[0] = ioe + 2*s_e[0] + sim.map(0, 0);  // e_0_0
            s_1[1] = ios + 3*i + 1;
            s_1[2] = ioe + 2*s_e[2] + sim.map(2, 1);  // e_2_1

            s_2[0] = ioe + 2*s_e[0] + sim.map(0, 1);  // e_0_1
            s_2[1] = ioe + 2*s_e[1] + sim.map(1, 0);  // e_1_0
            s_2[2] = ios + 3*i + 2;

            s_3[0] = ios + 3*i + 0;
            s_3[1] = ios + 3*i + 1;
            s_3[2] = ios + 3*i + 2;
          }

          // return fine triangle count
          return 4*num_simps;
        }
      }; // StandardIndexRefiner<Simplex<2>,2,1>

      /* ************************************************************************************* */
      /* ************************************************************************************* */
      /* ************************************************************************************* */

      /**
       * \brief StandardIndexRefiner implementation for Simplex<3>: Vertices-at-Edges
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Simplex<3>, 1, 0>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 1; // second template parameter
        static constexpr int face_dim = 0;  // third template parameter

        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index ioe = index_offsets[1];
          const Index ios = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set type (edges at simplices)
          typedef IndexSet<6> IndexSetTypeSE;

          // typedef index vector type
          typedef IndexSetTypeSE::IndexTupleType IndexTupleTypeSE;

          // fetch the edge-at-simplex index set
          const IndexSetTypeSE& index_set_s_e = index_set_holder_in.get_index_set<3,1>();

          // fetch number of coarse mesh simplices
          const Index num_simps = index_set_s_e.get_num_entities();

          /*                        v_3
                                    /\\
                                    /  \  \
                                  /    \    \
                                  /      \      \
                                /        \        \
                                /          \          \
                              /            \            \
                              /              \              x e_5
                            /                \                \
                            /                  \                  \
                          /                    \                    \
                          /                      \                      \ v_2
                    e_2 x                        x e_4             /    /
                        /                   s x    \           /        /
                      /                            \      /           /
                      /                              \ /              /
                    /                            /   \              /
                    /                        x         \            /
                  /                    /    e_1        \          x e_3
                  /                /                     \        /
                /            /                           \      /
                /        /                                 \    /
              /    /                                       \  /
              //_______________________x_____________________\/
            v_0                      e_0                    v_1

            orientation:
              t0: v_1 - v_2 - v_3
              t1: v_0 - v_2 - v_3
              t2: v_0 - v_1 - v_3
              t3: v_0 - v_1 - v_2
          */


          // loop over all coarse mesh simplices
          for(Index i(0); i < num_simps; ++i)
          {
            // fetch coarse mesh simplex-at-edge index vector
            const IndexTupleTypeSE& s_e = index_set_s_e[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexTupleType& e_0 = index_set_out[offset + 6*i + 0];
            IndexTupleType& e_1 = index_set_out[offset + 6*i + 1];
            IndexTupleType& e_2 = index_set_out[offset + 6*i + 2];
            IndexTupleType& e_3 = index_set_out[offset + 6*i + 3];
            IndexTupleType& e_4 = index_set_out[offset + 6*i + 4];
            IndexTupleType& e_5 = index_set_out[offset + 6*i + 5];

            e_0[0] = ioe + s_e[0];  //e_0
            e_0[1] = ios + i;       //s

            e_1[0] = ioe + s_e[1];  //e_1
            e_1[1] = ios + i;       //s

            e_2[0] = ioe + s_e[2];  //e_2
            e_2[1] = ios + i;       //s

            e_3[0] = ioe + s_e[3];  //e_3
            e_3[1] = ios + i;       //s

            e_4[0] = ioe + s_e[4];  //e_4
            e_4[1] = ios + i;       //s

            e_5[0] = ioe + s_e[5];  //e_5
            e_5[1] = ios + i;       //s
          }

          // return fine edge count
          return 6*num_simps;
        }
      }; // StandardIndexRefiner<Simplex<3>,1,0>

      /**
       * \brief StandardIndexRefiner implementation for Simplex<3>: Vertices-at-Tris
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Simplex<3>, 2, 0>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 2; // second template parameter
        static constexpr int face_dim = 0; // third template parameter

        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index ioe = index_offsets[1];
          const Index ios = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types
          typedef IndexSet<4> IndexSetTypeST;
          typedef IndexSet<6> IndexSetTypeSE;

          // typedef index vector type
          typedef IndexSetTypeSE::IndexTupleType IndexTupleTypeSE;

          // fetch the index sets
          const IndexSetTypeST& index_set_s_t = index_set_holder_in.get_index_set<3,2>();
          const IndexSetTypeSE& index_set_s_e = index_set_holder_in.get_index_set<3,1>();

          // fetch number of coarse mesh simplices
          const Index num_simps = index_set_s_t.get_num_entities();

          /*                        v_3
                                    /\\
                                    /  \  \
                                  /    \    \
                                  /      \      \
                                /        \        \
                                /          \          \
                              /            \           x e_5
                              /              \              \
                            /                \                \
                            /                  \                  \
                          /                    \                    \
                          /                  e_4 \                      \ v_2
                    e_2 x\                       x                 /    /
                        / \   \             s x    \           /        /
                      /   \      \                 \      /           /
                      /     \         ^              \ /              /
                    /       \            \  e_1  /   \              /
                    /         v    s_0        x        \            /
                  /           \        /    /          \          x e_3
                  /             \  /       /             \        /
                /            /  \       ^                \      /
                /        /        \    /                   \    /
              /    /              \ /                      \  /
              //____________________x________________________\/
            v_0                    e_0                      v_1

            orientation:
              t0: v_1 - v_2 - v_3
              t1: v_0 - v_2 - v_3
              t2: v_0 - v_1 - v_3
              t3: v_0 - v_1 - v_2
          */

          // loop over all coarse mesh simplices
          for(Index i(0); i < num_simps; ++i)
          {
            // fetch coarse mesh edges-at-simplex index vectors
            const IndexTupleTypeSE& s_e = index_set_s_e[i];

            // fetch fine mesh vertices-at-triangle index vectors
            IndexTupleType& t_0 = index_set_out[offset + 16*i + 0];
            IndexTupleType& t_1 = index_set_out[offset + 16*i + 1];
            IndexTupleType& t_2 = index_set_out[offset + 16*i + 2];
            IndexTupleType& t_3 = index_set_out[offset + 16*i + 3];
            IndexTupleType& t_4 = index_set_out[offset + 16*i + 4];
            IndexTupleType& t_5 = index_set_out[offset + 16*i + 5];
            IndexTupleType& t_6 = index_set_out[offset + 16*i + 6];
            IndexTupleType& t_7 = index_set_out[offset + 16*i + 7];
            IndexTupleType& t_8 = index_set_out[offset + 16*i + 8];
            IndexTupleType& t_9 = index_set_out[offset + 16*i + 9];
            IndexTupleType& t_10 = index_set_out[offset + 16*i + 10];
            IndexTupleType& t_11 = index_set_out[offset + 16*i + 11];
            IndexTupleType& t_12 = index_set_out[offset + 16*i + 12];
            IndexTupleType& t_13 = index_set_out[offset + 16*i + 13];
            IndexTupleType& t_14 = index_set_out[offset + 16*i + 14];
            IndexTupleType& t_15 = index_set_out[offset + 16*i + 15];

            // calculate fine simplex-vertex indices

            // triangles at the corners (e.g. s_0, the orientation is the same as for triangle(local point_index))
            t_0[0] = ioe + s_e[0]; // e_0
            t_0[1] = ioe + s_e[1]; // e_1
            t_0[2] = ioe + s_e[2]; // e_2

            t_1[0] = ioe + s_e[0]; // e_0
            t_1[1] = ioe + s_e[3]; // e_3
            t_1[2] = ioe + s_e[4]; // e_4

            t_2[0] = ioe + s_e[1]; // e_1
            t_2[1] = ioe + s_e[3]; // e_3
            t_2[2] = ioe + s_e[5]; // e_5

            t_3[0] = ioe + s_e[2]; // e_0
            t_3[1] = ioe + s_e[4]; // e_1
            t_3[2] = ioe + s_e[5]; // e_2

            // triangles adjuncted to the center (orientation i < j -> e_i - e_j - s)
            t_4[0] = ioe + s_e[0]; // e_0
            t_4[1] = ioe + s_e[1]; // e_1
            t_4[2] = ios + i;      // s

            t_5[0] = ioe + s_e[0]; // e_0
            t_5[1] = ioe + s_e[2]; // e_2
            t_5[2] = ios + i;      // s

            t_6[0] = ioe + s_e[0]; // e_0
            t_6[1] = ioe + s_e[3]; // e_3
            t_6[2] = ios + i;      // s

            t_7[0] = ioe + s_e[0]; // e_0
            t_7[1] = ioe + s_e[4]; // e_4
            t_7[2] = ios + i;      // s

            t_8[0] = ioe + s_e[1]; // e_1
            t_8[1] = ioe + s_e[2]; // e_2
            t_8[2] = ios + i;      // s

            t_9[0] = ioe + s_e[1]; // e_1
            t_9[1] = ioe + s_e[3]; // e_3
            t_9[2] = ios + i;      // s

            t_10[0] = ioe + s_e[1]; // e_1
            t_10[1] = ioe + s_e[5]; // e_5
            t_10[2] = ios + i;      // s

            t_11[0] = ioe + s_e[2]; // e_2
            t_11[1] = ioe + s_e[4]; // e_4
            t_11[2] = ios + i;      // s

            t_12[0] = ioe + s_e[2]; // e_2
            t_12[1] = ioe + s_e[5]; // e_5
            t_12[2] = ios + i;      // s

            t_13[0] = ioe + s_e[3]; // e_3
            t_13[1] = ioe + s_e[4]; // e_4
            t_13[2] = ios + i;      // s

            t_14[0] = ioe + s_e[3]; // e_3
            t_14[1] = ioe + s_e[5]; // e_5
            t_14[2] = ios + i;      // s

            t_15[0] = ioe + s_e[4]; // e_4
            t_15[1] = ioe + s_e[5]; // e_5
            t_15[2] = ios + i;      // s
          }

          // return fine triangle count
          return 16*num_simps;
        }
      }; // StandardIndexRefiner<Simplex<3>,2,0>

      /**
       * \brief StandardIndexRefiner implementation for Simplex<3>: Vertices-at-Tetras
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Simplex<3>, 3, 0>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 3; // second template parameter
        static constexpr int face_dim = 0;  // third template parameter

        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioe = index_offsets[1];
          const Index ios = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types
          typedef IndexSet<4> IndexSetTypeST;
          typedef IndexSet<6> IndexSetTypeSE;
          typedef IndexSet<4> IndexSetTypeSV;

          // typedef index vector type
          typedef IndexSetTypeSE::IndexTupleType IndexTupleTypeSE;
          typedef IndexSetTypeSV::IndexTupleType IndexTupleTypeSV;

          // fetch the index sets
          const IndexSetTypeST& index_set_s_t = index_set_holder_in.get_index_set<3,2>();
          const IndexSetTypeSE& index_set_s_e = index_set_holder_in.get_index_set<3,1>();
          const IndexSetTypeSV& index_set_s_v = index_set_holder_in.get_index_set<3,0>();

          // fetch number of coarse mesh simplices
          const Index num_simps = index_set_s_t.get_num_entities();

          /*                        v_3
                                    /\\
                                    /  \  \
                                  /    \    \
                                  /      \      \
                                /        \        \
                                /          \          \
                              /            \           x e_5
                              /              \              \
                            /                \                \
                            /                  \                  \
                          /                    \                    \
                          /                  e_4 \                      \ v_2
                    e_2 x\                       x                 /    /
                        / \   \             s x    \           /        /
                      /   \      \                 \      /           /
                      /     \         ^              \ /              /
                    /       \            \  e_1  /   \              /
                    /         v    s_0        x        \            /
                  /           \        /    /          \          x e_3
                  /             \  /       /             \        /
                /            /  \       ^                \      /
                /        /        \    /                   \    /
              /    /              \ /                      \  /
              //____________________x________________________\/
            v_0                    e_0                      v_1

            orientation:
              t0: v_1 - v_2 - v_3
              t1: v_0 - v_2 - v_3
              t2: v_0 - v_1 - v_3
              t3: v_0 - v_1 - v_2
          */

          // loop over all coarse mesh simplices
          for(Index i(0); i < num_simps; ++i)
          {
            // fetch coarse mesh index vectors
            const IndexTupleTypeSE& s_e = index_set_s_e[i];
            const IndexTupleTypeSV& s_v = index_set_s_v[i];

            // fetch fine mesh vertices-at-simplex index vectors
            IndexTupleType& s_0 = index_set_out[offset + 12*i + 0];
            IndexTupleType& s_1 = index_set_out[offset + 12*i + 1];
            IndexTupleType& s_2 = index_set_out[offset + 12*i + 2];
            IndexTupleType& s_3 = index_set_out[offset + 12*i + 3];
            IndexTupleType& s_4 = index_set_out[offset + 12*i + 4];
            IndexTupleType& s_5 = index_set_out[offset + 12*i + 5];
            IndexTupleType& s_6 = index_set_out[offset + 12*i + 6];
            IndexTupleType& s_7 = index_set_out[offset + 12*i + 7];
            IndexTupleType& s_8 = index_set_out[offset + 12*i + 8];
            IndexTupleType& s_9 = index_set_out[offset + 12*i + 9];
            IndexTupleType& s_10 = index_set_out[offset + 12*i + 10];
            IndexTupleType& s_11 = index_set_out[offset + 12*i + 11];

            // calculate fine simplex-vertex indices

            // simplices at the corners

            // v_0
            s_0[0] = ioe + s_e[0]; // e_0
            s_0[1] = ioe + s_e[2]; // e_2
            s_0[2] = ioe + s_e[1]; // e_1
            s_0[3] = iov + s_v[0]; // v_0

            s_1[0] = ioe + s_e[0]; // e_0
            s_1[1] = ioe + s_e[1]; // e_1
            s_1[2] = ioe + s_e[2]; // e_2
            s_1[3] = ios + i;      // s

            // v_1
            s_2[0] = ioe + s_e[0]; // e_0
            s_2[1] = ioe + s_e[3]; // e_3
            s_2[2] = ioe + s_e[4]; // e_4
            s_2[3] = iov + s_v[1]; // v_0

            s_3[0] = ioe + s_e[0]; // e_0
            s_3[1] = ioe + s_e[4]; // e_4
            s_3[2] = ioe + s_e[3]; // e_3
            s_3[3] = ios + i;      // s

            // v_2
            s_4[0] = ioe + s_e[1]; // e_1
            s_4[1] = ioe + s_e[5]; // e_5
            s_4[2] = ioe + s_e[3]; // e_3
            s_4[3] = iov + s_v[2]; // v_2

            s_5[0] = ioe + s_e[1]; // e_1
            s_5[1] = ioe + s_e[3]; // e_3
            s_5[2] = ioe + s_e[5]; // e_5
            s_5[3] = ios + i;      // s

            // v_3
            s_6[0] = ioe + s_e[2]; // e_2
            s_6[1] = ioe + s_e[4]; // e_4
            s_6[2] = ioe + s_e[5]; // e_5
            s_6[3] = iov + s_v[3]; // v_3

            s_7[0] = ioe + s_e[2]; // e_2
            s_7[1] = ioe + s_e[5]; // e_5
            s_7[2] = ioe + s_e[4]; // e_4
            s_7[3] = ios + i;      // s

            // simplices at the coarse mesh faces

            s_8[0] = ioe + s_e[3]; // e_3
            s_8[1] = ioe + s_e[4]; // e_4
            s_8[2] = ioe + s_e[5]; // e_5
            s_8[3] = ios + i;      // s

            s_9[0] = ioe + s_e[1]; // e_1
            s_9[1] = ioe + s_e[5]; // e_5
            s_9[2] = ioe + s_e[2]; // e_2
            s_9[3] = ios + i;      // s

            s_10[0] = ioe + s_e[0]; // e_0
            s_10[1] = ioe + s_e[2]; // e_2
            s_10[2] = ioe + s_e[4]; // e_4
            s_10[3] = ios + i;      // s

            s_11[0] = ioe + s_e[0]; // e_0
            s_11[1] = ioe + s_e[3]; // e_3
            s_11[2] = ioe + s_e[1]; // e_1
            s_11[3] = ios + i;      // s
          }
          // return fine simplex count
          return 12*num_simps;
        }
      }; // StandardIndexRefiner<Simplex<3>,3,0>

      /**
       * \brief StandardIndexRefiner implementation for Simplex<3>: Edges-at-Tris
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Simplex<3>, 2, 1>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 2; // second template parameter
        static constexpr int face_dim = 1; // third template parameter

        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch index offsets
          const Index iot = index_offsets[2];
          const Index ios = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types;
          typedef IndexSet<3> IndexSetTypeTV;
          typedef IndexSet<4> IndexSetTypeSV;
          typedef IndexSet<4> IndexSetTypeST;

          // typedef index vector type
          typedef IndexSetTypeSV::IndexTupleType IndexTupleTypeSV;
          typedef IndexSetTypeST::IndexTupleType IndexTupleTypeST;

          // fetch the index sets
          const IndexSetTypeSV& index_set_s_v = index_set_holder_in.get_index_set<3,0>();
          const IndexSetTypeST& index_set_s_t = index_set_holder_in.get_index_set<3,2>();
          const IndexSetTypeTV& index_set_t_v = index_set_holder_in.get_index_set<2,0>();

          // fetch number of coarse mesh simplices
          const Index num_simps = index_set_s_t.get_num_entities();

          // typedef the sub-index mapping
          typedef Intern::SubIndexMapping<ShapeType, cell_dim, face_dim> SubIndexMappingType;

          // loop over all coarse mesh simplices
          for(Index i(0); i < num_simps; ++i)
          {
            // fetch coarse mesh vertices-at-simplex and edges-at-simplex index vectors
            const IndexTupleTypeSV& s_v = index_set_s_v[i];
            const IndexTupleTypeST& s_t = index_set_s_t[i];

            // create a sub-index mapping object
            SubIndexMappingType sim(s_v, s_t, index_set_t_v);

            // fetch fine mesh edges-at-triangle index vectors
            IndexTupleType& t_0 = index_set_out[offset + 16*i + 0];
            IndexTupleType& t_1 = index_set_out[offset + 16*i + 1];
            IndexTupleType& t_2 = index_set_out[offset + 16*i + 2];
            IndexTupleType& t_3 = index_set_out[offset + 16*i + 3];
            IndexTupleType& t_4 = index_set_out[offset + 16*i + 4];
            IndexTupleType& t_5 = index_set_out[offset + 16*i + 5];
            IndexTupleType& t_6 = index_set_out[offset + 16*i + 6];
            IndexTupleType& t_7 = index_set_out[offset + 16*i + 7];
            IndexTupleType& t_8 = index_set_out[offset + 16*i + 8];
            IndexTupleType& t_9 = index_set_out[offset + 16*i + 9];
            IndexTupleType& t_10 = index_set_out[offset + 16*i + 10];
            IndexTupleType& t_11 = index_set_out[offset + 16*i + 11];
            IndexTupleType& t_12 = index_set_out[offset + 16*i + 12];
            IndexTupleType& t_13 = index_set_out[offset + 16*i + 13];
            IndexTupleType& t_14 = index_set_out[offset + 16*i + 14];
            IndexTupleType& t_15 = index_set_out[offset + 16*i + 15];

            // calculate fine edges-at-triangle indices

            // triangle at the corners
            t_0[0] = iot + 3*s_t[1] + sim.map(1, 0);
            t_0[1] = iot + 3*s_t[2] + sim.map(2, 0);
            t_0[2] = iot + 3*s_t[3] + sim.map(3, 0);

            t_1[0] = iot + 3*s_t[0] + sim.map(0, 0);
            t_1[1] = iot + 3*s_t[2] + sim.map(2, 1);
            t_1[2] = iot + 3*s_t[3] + sim.map(3, 1);

            t_2[0] = iot + 3*s_t[0] + sim.map(0, 1);
            t_2[1] = iot + 3*s_t[1] + sim.map(1, 1);
            t_2[2] = iot + 3*s_t[3] + sim.map(3, 2);

            t_3[0] = iot + 3*s_t[0] + sim.map(0, 2);
            t_3[1] = iot + 3*s_t[1] + sim.map(1, 2);
            t_3[2] = iot + 3*s_t[2] + sim.map(2, 2);

            // triangles at the midpoint

            //e_0 triangles
            t_4[0] = ios + 6*i + 1;
            t_4[1] = ios + 6*i + 0;
            t_4[2] = iot + 3*s_t[3] + sim.map(3, 0);

            t_5[0] = ios + 6*i + 2;
            t_5[1] = ios + 6*i + 0;
            t_5[2] = iot + 3*s_t[2] + sim.map(2, 0);

            t_6[0] = ios + 6*i + 3;
            t_6[1] = ios + 6*i + 0;
            t_6[2] = iot + 3*s_t[3] + sim.map(3, 1);

            t_7[0] = ios + 6*i + 4;
            t_7[1] = ios + 6*i + 0;
            t_7[2] = iot + 3*s_t[2] + sim.map(2, 1);

            //e_1 triangles
            t_8[0] = ios + 6*i + 2;
            t_8[1] = ios + 6*i + 1;
            t_8[2] = iot + 3*s_t[1] + sim.map(1, 0);

            t_9[0] = ios + 6*i + 3;
            t_9[1] = ios + 6*i + 1;
            t_9[2] = iot + 3*s_t[3] + sim.map(3, 2);

            t_10[0] = ios + 6*i + 5;
            t_10[1] = ios + 6*i + 1;
            t_10[2] = iot + 3*s_t[1] + sim.map(1, 1);

            //e_2 triangles
            t_11[0] = ios + 6*i + 4;
            t_11[1] = ios + 6*i + 2;
            t_11[2] = iot + 3*s_t[2] + sim.map(2, 2);

            t_12[0] = ios + 6*i + 5;
            t_12[1] = ios + 6*i + 2;
            t_12[2] = iot + 3*s_t[1] + sim.map(1, 2);

            //e_3 triangles
            t_13[0] = ios + 6*i + 4;
            t_13[1] = ios + 6*i + 3;
            t_13[2] = iot + 3*s_t[0] + sim.map(0, 0);

            t_14[0] = ios + 6*i + 5;
            t_14[1] = ios + 6*i + 3;
            t_14[2] = iot + 3*s_t[0] + sim.map(0, 1);

            //e_4 triangles
            t_15[0] = ios + 6*i + 5;
            t_15[1] = ios + 6*i + 4;
            t_15[2] = iot + 3*s_t[0] + sim.map(0, 2);

          }
          // return fine triangle count
          return 16*num_simps;
        }
      }; // StandardIndexRefiner<Simplex<3>,2,1>

      /**
       * \brief StandardIndexRefiner implementation for Simplex<3>: Edges-at-Tetras
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Simplex<3>, 3, 1>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 3; // second template parameter
        static constexpr int face_dim = 1; // third template parameter

        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch index offsets
          const Index ioe = index_offsets[1];
          const Index iot = index_offsets[2];
          const Index ios = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types
          typedef IndexSet<2> IndexSetTypeEV;
          typedef IndexSet<3> IndexSetTypeTV;
          typedef IndexSet<4> IndexSetTypeSV;
          typedef IndexSet<6> IndexSetTypeSE;
          typedef IndexSet<4> IndexSetTypeST;

          // typedef index vector type
          typedef IndexSetTypeSV::IndexTupleType IndexTupleTypeSV;
          typedef IndexSetTypeSE::IndexTupleType IndexTupleTypeSE;
          typedef IndexSetTypeST::IndexTupleType IndexTupleTypeST;

          // fetch the index sets
          const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();
          const IndexSetTypeTV& index_set_t_v = index_set_holder_in.get_index_set<2,0>();
          const IndexSetTypeSV& index_set_s_v = index_set_holder_in.get_index_set<3,0>();
          const IndexSetTypeSE& index_set_s_e = index_set_holder_in.get_index_set<3,1>();
          const IndexSetTypeST& index_set_s_t = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh simplices
          const Index num_simps = index_set_s_t.get_num_entities();

          // typedef the sub-index mapping
          typedef Intern::SubIndexMapping<ShapeType, 2, 1> SubIndexMappingType;
          typedef Intern::SubIndexMapping<ShapeType, face_dim, 0> EdgeIndexMappingType;

          // loop over all coarse mesh simplices
          for(Index i(0); i < num_simps; ++i)
          {
            // fetch coarse mesh index vectors
            const IndexTupleTypeSV& s_v = index_set_s_v[i];
            const IndexTupleTypeSE& s_e = index_set_s_e[i];
            const IndexTupleTypeST& s_t = index_set_s_t[i];

            // create sub-index mapping objects
            SubIndexMappingType sim(s_v, s_t, index_set_t_v);
            EdgeIndexMappingType edgesim(s_v, s_e, index_set_e_v);

            // fetch fine mesh edges-at-simplex index vectors
            IndexTupleType& s_0 = index_set_out[offset + 12*i + 0];
            IndexTupleType& s_1 = index_set_out[offset + 12*i + 1];
            IndexTupleType& s_2 = index_set_out[offset + 12*i + 2];
            IndexTupleType& s_3 = index_set_out[offset + 12*i + 3];
            IndexTupleType& s_4 = index_set_out[offset + 12*i + 4];
            IndexTupleType& s_5 = index_set_out[offset + 12*i + 5];
            IndexTupleType& s_6 = index_set_out[offset + 12*i + 6];
            IndexTupleType& s_7 = index_set_out[offset + 12*i + 7];
            IndexTupleType& s_8 = index_set_out[offset + 12*i + 8];
            IndexTupleType& s_9 = index_set_out[offset + 12*i + 9];
            IndexTupleType& s_10 = index_set_out[offset + 12*i + 10];
            IndexTupleType& s_11 = index_set_out[offset + 12*i + 11];

            // calculate fine edges-at-simplex indices

            // v_0-simplices
            s_0[0] = iot + 3*s_t[2] + sim.map(2, 0);
            s_0[1] = iot + 3*s_t[3] + sim.map(3, 0);
            s_0[2] = ioe + 2*s_e[0] + edgesim.map(0, 0);
            s_0[3] = iot + 3*s_t[1] + sim.map(1, 0);
            s_0[4] = ioe + 2*s_e[2] + edgesim.map(2, 0);
            s_0[5] = ioe + 2*s_e[1] + edgesim.map(1, 0);

            s_1[0] = iot + 3*s_t[3] + sim.map(3, 0);
            s_1[1] = iot + 3*s_t[2] + sim.map(2, 0);
            s_1[2] = ios + 6*i + 0;
            s_1[3] = iot + 3*s_t[1] + sim.map(1, 0);
            s_1[4] = ios + 6*i + 1;
            s_1[5] = ios + 6*i + 2;

            // v_1-simplices
            s_2[0] = iot + 3*s_t[3] + sim.map(3, 1);
            s_2[1] = iot + 3*s_t[2] + sim.map(2, 1);
            s_2[2] = ioe + 2*s_e[0] + edgesim.map(0, 1);
            s_2[3] = iot + 3*s_t[0] + sim.map(0, 0);
            s_2[4] = ioe + 2*s_e[3] + edgesim.map(3, 0);
            s_2[5] = ioe + 2*s_e[4] + edgesim.map(4, 0);

            s_3[0] = iot + 3*s_t[2] + sim.map(2, 1);
            s_3[1] = iot + 3*s_t[3] + sim.map(3, 1);
            s_3[2] = ios + 6*i + 0;
            s_3[3] = iot + 3*s_t[0] + sim.map(0, 0);
            s_3[5] = ios + 6*i + 3;
            s_3[4] = ios + 6*i + 4;

            // v_2-simplices
            s_4[0] = iot + 3*s_t[1] + sim.map(1, 1);
            s_4[1] = iot + 3*s_t[3] + sim.map(3, 2);
            s_4[2] = ioe + 2*s_e[1] + edgesim.map(1, 1);
            s_4[3] = iot + 3*s_t[0] + sim.map(0, 1);
            s_4[4] = ioe + 2*s_e[5] + edgesim.map(5, 0);
            s_4[5] = ioe + 2*s_e[3] + edgesim.map(3, 1);

            s_5[0] = iot + 3*s_t[3] + sim.map(3, 2);
            s_5[1] = iot + 3*s_t[1] + sim.map(1, 1);
            s_5[2] = ios + 6*i + 1;
            s_5[3] = iot + 3*s_t[0] + sim.map(0, 1);
            s_5[4] = ios + 6*i + 3;
            s_5[5] = ios + 6*i + 5;

            // v_3-simplices
            s_6[0] = iot + 3*s_t[2] + sim.map(2, 2);
            s_6[1] = iot + 3*s_t[1] + sim.map(1, 2);
            s_6[2] = ioe + 2*s_e[2] + edgesim.map(2, 1);
            s_6[3] = iot + 3*s_t[0] + sim.map(0, 2);
            s_6[4] = ioe + 2*s_e[4] + edgesim.map(4, 1);
            s_6[5] = ioe + 2*s_e[5] + edgesim.map(5, 1);

            s_7[0] = iot + 3*s_t[1] + sim.map(1, 2);
            s_7[1] = iot + 3*s_t[2] + sim.map(2, 2);
            s_7[2] = ios + 6*i + 2;
            s_7[3] = iot + 3*s_t[0] + sim.map(0, 2);
            s_7[4] = ios + 6*i + 5;
            s_7[5] = ios + 6*i + 4;

            //triangle simplices

            s_8[0] = iot + 3*s_t[0] + sim.map(0, 0);
            s_8[1] = iot + 3*s_t[0] + sim.map(0, 1);
            s_8[2] = ios + 6*i + 3;
            s_8[3] = iot + 3*s_t[0] + sim.map(0, 2);
            s_8[4] = ios + 6*i + 4;
            s_8[5] = ios + 6*i + 5;

            s_9[0] = iot + 3*s_t[1] + sim.map(1, 1);
            s_9[1] = iot + 3*s_t[1] + sim.map(1, 0);
            s_9[2] = ios + 6*i + 1;
            s_9[3] = iot + 3*s_t[1] + sim.map(1, 2);
            s_9[4] = ios + 6*i + 5;
            s_9[5] = ios + 6*i + 2;

            s_10[0] = iot + 3*s_t[2] + sim.map(2, 0);
            s_10[1] = iot + 3*s_t[2] + sim.map(2, 1);
            s_10[2] = ios + 6*i + 0;
            s_10[3] = iot + 3*s_t[2] + sim.map(2, 2);
            s_10[4] = ios + 6*i + 2;
            s_10[5] = ios + 6*i + 4;

            s_11[0] = iot + 3*s_t[3] + sim.map(3, 1);
            s_11[1] = iot + 3*s_t[3] + sim.map(3, 0);
            s_11[2] = ios + 6*i + 0;
            s_11[3] = iot + 3*s_t[3] + sim.map(3, 2);
            s_11[4] = ios + 6*i + 3;
            s_11[5] = ios + 6*i + 1;

          }
          // return fine simplex count
          return 12*num_simps;
        }
      }; // StandardIndexRefiner<Simplex<3>,3,1>

      /**
       * \brief StandardIndexRefiner implementation for Simplex<3>: Tris-at-Tetras
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Simplex<3>, 3, 2>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 3; // second template parameter
        static constexpr int face_dim = 2; // third template parameter

        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch index offsets
          const Index iot = index_offsets[2];
          const Index ios = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types
          typedef IndexSet<3> IndexSetTypeTV;
          typedef IndexSet<4> IndexSetTypeSV;
          typedef IndexSet<4> IndexSetTypeST;

          // typedef index vector type
          typedef IndexSetTypeSV::IndexTupleType IndexTupleTypeSV;
          typedef IndexSetTypeST::IndexTupleType IndexTupleTypeST;

          // fetch the index sets
          const IndexSetTypeTV& index_set_t_v = index_set_holder_in.get_index_set<2,0>();
          const IndexSetTypeSV& index_set_s_v = index_set_holder_in.get_index_set<3,0>();
          const IndexSetTypeST& index_set_s_t = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh simplices
          const Index num_simps = index_set_s_t.get_num_entities();

          // typedef the sub-index mapping
          typedef Intern::SubIndexMapping<ShapeType, face_dim, 0> SubIndexMappingType;

          // loop over all coarse mesh simplices
          for(Index i(0); i < num_simps; ++i)
          {
            // fetch coarse mesh vertices-at-simplex and triangle-at-simplex index vectors
            const IndexTupleTypeSV& s_v = index_set_s_v[i];
            const IndexTupleTypeST& s_t = index_set_s_t[i];

            // create a sub-index mapping object
            SubIndexMappingType sim(s_v, s_t, index_set_t_v);

            // fetch fine mesh triangle-at-simplex index vectors
            IndexTupleType& s_0 = index_set_out[offset + 12*i + 0];
            IndexTupleType& s_1 = index_set_out[offset + 12*i + 1];
            IndexTupleType& s_2 = index_set_out[offset + 12*i + 2];
            IndexTupleType& s_3 = index_set_out[offset + 12*i + 3];
            IndexTupleType& s_4 = index_set_out[offset + 12*i + 4];
            IndexTupleType& s_5 = index_set_out[offset + 12*i + 5];
            IndexTupleType& s_6 = index_set_out[offset + 12*i + 6];
            IndexTupleType& s_7 = index_set_out[offset + 12*i + 7];
            IndexTupleType& s_8 = index_set_out[offset + 12*i + 8];
            IndexTupleType& s_9 = index_set_out[offset + 12*i + 9];
            IndexTupleType& s_10 = index_set_out[offset + 12*i + 10];
            IndexTupleType& s_11 = index_set_out[offset + 12*i + 11];

            // calculate fine triangle-at-simplex indices

            // v_0-cube
            s_0[0] = iot + 4*s_t[1] + sim.map(1, 0);
            s_0[2] = iot + 4*s_t[2] + sim.map(2, 0);
            s_0[1] = iot + 4*s_t[3] + sim.map(3, 0);
            s_0[3] = ios + 16*i + 0;

            s_1[0] = ios + 16*i + 8;
            s_1[1] = ios + 16*i + 5;
            s_1[2] = ios + 16*i + 4;
            s_1[3] = ios + 16*i + 0;

            // v_1-cube
            s_2[0] = iot + 4*s_t[0] + sim.map(0, 0);
            s_2[1] = iot + 4*s_t[2] + sim.map(2, 1);
            s_2[2] = iot + 4*s_t[3] + sim.map(3, 1);
            s_2[3] = ios + 16*i + 1;

            s_3[0] = ios + 16*i + 13;
            s_3[2] = ios + 16*i + 7;
            s_3[1] = ios + 16*i + 6;
            s_3[3] = ios + 16*i + 1;

            // v_2-cube
            s_4[0] = iot + 4*s_t[0] + sim.map(0, 1);
            s_4[2] = iot + 4*s_t[1] + sim.map(1, 1);
            s_4[1] = iot + 4*s_t[3] + sim.map(3, 2);
            s_4[3] = ios + 16*i + 2;

            s_5[0] = ios + 16*i + 14;
            s_5[1] = ios + 16*i + 10;
            s_5[2] = ios + 16*i + 9;
            s_5[3] = ios + 16*i + 2;

            // v_3-cube
            s_6[0] = iot + 4*s_t[0] + sim.map(0, 2);
            s_6[1] = iot + 4*s_t[1] + sim.map(1, 2);
            s_6[2] = iot + 4*s_t[2] + sim.map(2, 2);
            s_6[3] = ios + 16*i + 3;

            s_7[0] = ios + 16*i + 15;
            s_7[2] = ios + 16*i + 12;
            s_7[1] = ios + 16*i + 11;
            s_7[3] = ios + 16*i + 3;

            //triangle simplices

            s_8[0] = ios + 16*i + 15;
            s_8[1] = ios + 16*i + 14;
            s_8[2] = ios + 16*i + 13;
            s_8[3] = iot + 4*s_t[0] + 3;

            s_9[0] = ios + 16*i + 12;
            s_9[1] = ios + 16*i + 8;
            s_9[2] = ios + 16*i + 10;
            s_9[3] = iot + 4*s_t[1] + 3;

            s_10[0] = ios + 16*i + 11;
            s_10[1] = ios + 16*i + 7;
            s_10[2] = ios + 16*i + 5;
            s_10[3] = iot + 4*s_t[2] + 3;

            s_11[0] = ios + 16*i + 9;
            s_11[1] = ios + 16*i + 4;
            s_11[2] = ios + 16*i + 6;
            s_11[3] = iot + 4*s_t[3] + 3;

          }
          // return fine simplex count
          return 12*num_simps;
        }
      }; // StandardIndexRefiner<Simplex<3>,3,2>
      /**
        * \brief StandardIndexRefiner implementation for Hypercube<1>: Vertices-at-Edges
        *
        * \author Peter Zajac
        */
      template<>
      struct StandardIndexRefiner<Shape::Hypercube<1>, 1, 0>
      {
        static constexpr int shape_dim = 1;
        static constexpr int cell_dim = 1; // second template parameter
        static constexpr int face_dim = 0; // third template parameter

        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioe = index_offsets[1];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set type; an edge has 2 vertices, so we need an IndexSet<2>
          typedef IndexSet<2> IndexSetTypeEV;

          // typedef index vector reference type
          typedef IndexSetTypeEV::IndexTupleType IndexTupleTypeEV;

          // fetch the vertices-at-edge index set
          const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();

          // fetch number of coarse mesh edges
          const Index num_edges = index_set_e_v.get_num_entities();

          // loop over all coarse mesh edges
          for(Index i(0); i < num_edges; ++i)
          {
            // fetch coarse mesh vertices-at-edge index vector
            const IndexTupleTypeEV& e_v = index_set_e_v[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexTupleType& idx0 = index_set_out[offset + 2*i + 0];
            IndexTupleType& idx1 = index_set_out[offset + 2*i + 1];

            // calculate fine edge-vertex indices
            idx0[0] = iov + e_v[0];
            idx0[1] = ioe + i;
            idx1[0] = ioe + i;
            idx1[1] = iov + e_v[1];
          }

          // return fine edge count
          return 2*num_edges;
        }
      }; // StandardIndexRefiner<Hypercube<1>,1,0>

      /* ************************************************************************************* */
      /* ************************************************************************************* */
      /* ************************************************************************************* */

      /**
       * \brief StandardIndexRefiner implementation for Hypercube<2>: Vertices-at-Edges
       *
       * \author Peter Zajac
       */
      template<>
      struct StandardIndexRefiner<Shape::Hypercube<2>, 1, 0>
      {
        static constexpr int shape_dim = 2;
        static constexpr int cell_dim = 1; // second template parameter
        static constexpr int face_dim = 0; // third template parameter

        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index ioe = index_offsets[1];
          const Index ioq = index_offsets[2];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set type; a quad has 4 edges, so we need an IndexSet<4>
          typedef IndexSet<4> IndexSetTypeQE;

          // typedef index vector type
          typedef IndexSetTypeQE::IndexTupleType IndexTupleTypeQE;

          // fetch the edges-at-quad index set
          const IndexSetTypeQE& index_set_q_e = index_set_holder_in.get_index_set<2,1>();

          // fetch number of coarse mesh quads
          const Index num_quads = index_set_q_e.get_num_entities();

          // Each coarse mesh quad generates four fine mesh edges (e_i) upon refinement, with
          // each one connecting one of the fine mesh vertices generated from edge refinement
          // (v_i) with the vertex generated from quad refinement (x).
          //
          //  +------------v_1------------+
          //  |             ^             |
          //  |             |             |
          //  |            e_1            |
          //  |             |             |
          //  |             |             |
          // v_2----e_2---->x-----e_3--->v_3
          //  |             ^             |
          //  |             |             |
          //  |            e_0            |
          //  |             |             |
          //  |             |             |
          //  +------------v_0------------+


          // loop over all coarse mesh edges
          for(Index i(0); i < num_quads; ++i)
          {
            // fetch coarse mesh edges-at-quad index vector
            const IndexTupleTypeQE& q_e = index_set_q_e[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexTupleType& e_0 = index_set_out[offset + 4*i + 0];
            IndexTupleType& e_1 = index_set_out[offset + 4*i + 1];
            IndexTupleType& e_2 = index_set_out[offset + 4*i + 2];
            IndexTupleType& e_3 = index_set_out[offset + 4*i + 3];

            e_0[0] = ioe + q_e[0]; // v_0
            e_0[1] = ioq + i;      // x

            e_1[0] = ioq + i;      // x
            e_1[1] = ioe + q_e[1]; // v_1

            e_2[0] = ioe + q_e[2]; // v_2
            e_2[1] = ioq + i;      // x

            e_3[0] = ioq + i;      // x
            e_3[1] = ioe + q_e[3]; // v_3
          }

          // return fine edge count
          return 4*num_quads;
        }
      }; // StandardIndexRefiner<Hypercube<2>,1,0>

      /**
       * \brief StandardIndexRefiner implementation for Hypercube<2>: Vertices-at-Quads
       *
       * \author Peter Zajac
       */
      template<>
      struct StandardIndexRefiner<Shape::Hypercube<2>, 2, 0>
      {
        static constexpr int shape_dim = 2;
        static constexpr int cell_dim = 2; // second template parameter
        static constexpr int face_dim = 0; // third template parameter

        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioe = index_offsets[1];
          const Index ioq = index_offsets[2];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types; a quad has 4 vertices and 4 edges, so both need an IndexSet<4>
          typedef IndexSet<4> IndexSetTypeQV;
          typedef IndexSet<4> IndexSetTypeQE;

          // typedef index vector type
          typedef IndexSetTypeQV::IndexTupleType IndexTupleTypeQV;
          typedef IndexSetTypeQE::IndexTupleType IndexTupleTypeQE;

          // fetch the quad-vertex- and the quad-edge index set
          const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
          const IndexSetTypeQE& index_set_q_e = index_set_holder_in.get_index_set<2,1>();

          // fetch number of coarse mesh quads
          const Index num_quads = index_set_q_v.get_num_entities();

          // Each coarse mesh quad generates four fine mesh quads (q_i) upon refinement, with
          // each one connecting one of the coarse mesh vertices (v_i), its two neighbour fine
          // mesh vertices generated from edge refinement (e_i) and the vertex generated from
          // quad refinement (x):
          //
          // v_2-----------e_1-----------v_3
          //  |             ^             |
          //  |             |             |
          //  |     q_2     |     q_3     |
          //  |             |             |
          //  |             |             |
          // e_2----------->x----------->e_3
          //  |             ^             |
          //  |             |             |
          //  |     q_0     |     q_1     |
          //  |             |             |
          //  |             |             |
          // v_0-----------e_0-----------v_1

          // loop over all coarse mesh quads
          for(Index i(0); i < num_quads; ++i)
          {
            // fetch coarse mesh vertices-at-quad and edges-at-quad index vectors
            const IndexTupleTypeQV& q_v = index_set_q_v[i];
            const IndexTupleTypeQE& q_e = index_set_q_e[i];

            // fetch fine mesh vertices-at-quad index vectors
            IndexTupleType& q_0 = index_set_out[offset + 4*i + 0];
            IndexTupleType& q_1 = index_set_out[offset + 4*i + 1];
            IndexTupleType& q_2 = index_set_out[offset + 4*i + 2];
            IndexTupleType& q_3 = index_set_out[offset + 4*i + 3];

            // calculate fine quad-vertex indices
            q_0[0] = iov + q_v[0]; // v_0
            q_0[1] = ioe + q_e[0]; // e_0
            q_0[2] = ioe + q_e[2]; // e_2
            q_0[3] = ioq + i;      // x

            q_1[0] = ioe + q_e[0]; // e_0
            q_1[1] = iov + q_v[1]; // v_1
            q_1[2] = ioq + i;      // x
            q_1[3] = ioe + q_e[3]; // e_3

            q_2[0] = ioe + q_e[2]; // e_2
            q_2[1] = ioq + i;      // x
            q_2[2] = iov + q_v[2]; // v_2
            q_2[3] = ioe + q_e[1]; // e_1

            q_3[0] = ioq + i;      // x
            q_3[1] = ioe + q_e[3]; // e_3
            q_3[2] = ioe + q_e[1]; // e_1
            q_3[3] = iov + q_v[3]; // v_3
          }

          // return fine quad count
          return 4*num_quads;
        }
      }; // StandardIndexRefiner<Hypercube<2>,2,0>

      /**
       * \brief StandardIndexRefiner implementation for Hypercube<2>: Edges-at-Quads
       *
       * \author Peter Zajac
       */
      template<>
      struct StandardIndexRefiner<Shape::Hypercube<2>, 2, 1>
      {
        static constexpr int shape_dim = 2;
        static constexpr int cell_dim = 2; // second template parameter
        static constexpr int face_dim = 1; // third template parameter

        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch edge index offsets
          const Index ioe = index_offsets[1];
          const Index ioq = index_offsets[2];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types;
          typedef IndexSet<2> IndexSetTypeEV; // an edge has 2 vertices
          typedef IndexSet<4> IndexSetTypeQV; // a quad has 4 vertices
          typedef IndexSet<4> IndexSetTypeQE; // a quad has 4 edges

          // typedef index vector type
          typedef IndexSetTypeQV::IndexTupleType IndexTupleTypeQV;
          typedef IndexSetTypeQE::IndexTupleType IndexTupleTypeQE;

          // fetch the index sets
          const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();
          const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
          const IndexSetTypeQE& index_set_q_e = index_set_holder_in.get_index_set<2,1>();

          // fetch number of coarse mesh quads
          const Index num_quads = index_set_q_v.get_num_entities();

          // typedef the sub-index mapping
          typedef Intern::SubIndexMapping<ShapeType, face_dim, 0> SubIndexMappingType;

          /// \todo some text for fancy ASCII-art
          //
          //   +----e_1_0----+----e_1_1----+
          //   |             |             |
          //   |             |             |
          // e_2_1   q_2    f_1    q_3   e_3_1
          //   |             |             |
          //   |             |             |
          //   +-----f_2-----+-----f_3-----+
          //   |             |             |
          //   |             |             |
          // e_2_0   q_0    f_0    q_1   e_3_0
          //   |             |             |
          //   |             |             |
          //   +----e_0_0----+----e_0_1----+

          // loop over all coarse mesh quads
          for(Index i(0); i < num_quads; ++i)
          {
            // fetch coarse mesh vertices-at-quad and edges-at-quad index vectors
            const IndexTupleTypeQV& q_v = index_set_q_v[i];
            const IndexTupleTypeQE& q_e = index_set_q_e[i];

            // create a sub-index mapping object
            SubIndexMappingType sim(q_v, q_e, index_set_e_v);

            // fetch fine mesh edges-at-quad index vectors
            IndexTupleType& q_0 = index_set_out[offset + 4*i + 0];
            IndexTupleType& q_1 = index_set_out[offset + 4*i + 1];
            IndexTupleType& q_2 = index_set_out[offset + 4*i + 2];
            IndexTupleType& q_3 = index_set_out[offset + 4*i + 3];

            // calculate fine edges-at-quad indices
            q_0[0] = ioe + 2*q_e[0] + sim.map(0, 0);  // e_0_0
            q_0[1] = ioq + 4*i + 2;                   // f_2
            q_0[2] = ioe + 2*q_e[2] + sim.map(2, 0);  // e_2_0
            q_0[3] = ioq + 4*i + 0;                   // f_0

            q_1[0] = ioe + 2*q_e[0] + sim.map(0, 1);  // e_0_1
            q_1[1] = ioq + 4*i + 3;                   // f_3
            q_1[2] = ioq + 4*i + 0;                   // f_0
            q_1[3] = ioe + 2*q_e[3] + sim.map(3, 0);  // e_3_0

            q_2[0] = ioq + 4*i + 2;                   // f_2
            q_2[1] = ioe + 2*q_e[1] + sim.map(1, 0);  // e_1_0
            q_2[2] = ioe + 2*q_e[2] + sim.map(2, 1);  // e_2_1
            q_2[3] = ioq + 4*i + 1;                   // f_1

            q_3[0] = ioq + 4*i + 3;                   // f_3
            q_3[1] = ioe + 2*q_e[1] + sim.map(1, 1);  // e_1_1
            q_3[2] = ioq + 4*i + 1;                   // f_1
            q_3[3] = ioe + 2*q_e[3] + sim.map(3, 1);  // e_3_1
          }

          // return fine quad count
          return 4*num_quads;
        }
      }; // StandardIndexRefiner<Hypercube<2>,2,1>

      /* ************************************************************************************* */
      /* ************************************************************************************* */
      /* ************************************************************************************* */

      /**
       * \brief StandardIndexRefiner implementation for Hypercube<3>: Vertices-at-Edges
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Hypercube<3>, 1, 0>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 1; // second template parameter
        static constexpr int face_dim = 0; // third template parameter

        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index ioq = index_offsets[2];
          const Index ioc = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set type
          typedef IndexSet<6> IndexSetTypeCQ;

          // typedef index vector type
          typedef IndexSetTypeCQ::IndexTupleType IndexTupleTypeCQ;

          // fetch the quad-at-cube index set
          const IndexSetTypeCQ& index_set_c_q = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh cubes
          const Index num_cubes = index_set_c_q.get_num_entities();

          //
          //                q_1
          //                 ^
          //                 |   q_3
          //               e_1  /
          //                 | e_3
          //                 |/
          // q_4-----e_4---->x-----e_5---->q_5
          //                /^
          //             e_2 |
          //              /  e_0
          //           q_2   |
          //                 |
          //                q_0


          // loop over all coarse mesh cubes
          for(Index i(0); i < num_cubes; ++i)
          {
            // fetch coarse mesh quad-at-cube index vector
            const IndexTupleTypeCQ& c_q = index_set_c_q[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexTupleType& e_0 = index_set_out[offset + 6*i + 0];
            IndexTupleType& e_1 = index_set_out[offset + 6*i + 1];
            IndexTupleType& e_2 = index_set_out[offset + 6*i + 2];
            IndexTupleType& e_3 = index_set_out[offset + 6*i + 3];
            IndexTupleType& e_4 = index_set_out[offset + 6*i + 4];
            IndexTupleType& e_5 = index_set_out[offset + 6*i + 5];

            e_0[0] = ioq + c_q[0]; // q_0
            e_0[1] = ioc + i;      // x

            e_1[0] = ioc + i;      // x
            e_1[1] = ioq + c_q[1]; // q_1

            e_2[0] = ioq + c_q[2]; // q_2
            e_2[1] = ioc + i;      // x

            e_3[0] = ioc + i;      // x
            e_3[1] = ioq + c_q[3]; // q_3

            e_4[0] = ioq + c_q[4]; // q_4
            e_4[1] = ioc + i;      // x

            e_5[0] = ioc + i;      // x
            e_5[1] = ioq + c_q[5]; // q_5
          }

          // return fine edge count
          return 6*num_cubes;
        }
      }; // StandardIndexRefiner<Hypercube<3>,1,0>

      /**
       * \brief StandardIndexRefiner implementation for Hypercube<3>: Vertices-at-Quads
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Hypercube<3>, 2, 0>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 2; // second template parameter
        static constexpr int face_dim = 0; // third template parameter

        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index ioe = index_offsets[1];
          const Index ioq = index_offsets[2];
          const Index ioc = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types
          typedef IndexSet<6> IndexSetTypeCQ;
          typedef IndexSet<12> IndexSetTypeCE;

          // typedef index vector type
          typedef IndexSetTypeCQ::IndexTupleType IndexTupleTypeCQ;
          typedef IndexSetTypeCE::IndexTupleType IndexTupleTypeCE;

          // fetch the index sets
          const IndexSetTypeCQ& index_set_c_q = index_set_holder_in.get_index_set<3,2>();
          const IndexSetTypeCE& index_set_c_e = index_set_holder_in.get_index_set<3,1>();

          // fetch number of coarse mesh cubes
          const Index num_cubes = index_set_c_q.get_num_entities();

          // x-normal surfaces
          //                         e_3
          //                         /|
          //                        / |
          //                       /  |
          //                    q_1 s |
          //                     /| 3 |
          //                    / |   q_3
          //                   /  |  /|
          //                e_2 s | / |
          //                  | 2 |/  |
          //                  |   x s |
          //                  |  /| 1 |
          //                  | / |   e_1
          //                  |/  |  /
          //                q_2 s | /
          //                  | 0 |/
          //                  |  q_0
          //                  |  /                     z
          //                  | /                      ^ y
          //                  |/                       |/
          //                 e_0                       ---->x


          //  y-normal surfaces
          //
          //       e_6-----------q_1-----------e_7
          //        |             |             |
          //        |             |             |
          //        |     s_6     |     s_7     |
          //        |             |             |
          //        |             |             |
          //       q_4------------x------------q_5
          //        |             |             |
          //        |             |             |
          //        |     s_4     |     s_5     |      z
          //        |             |             |      ^ y
          //        |             |             |      |/
          //       e_4-----------q_0-----------e_5     ---->x


          // z-normal surfaces
          //
          //           e_10----------q_3----------e_11
          //            /            /            /
          //           /   s_10     /   s_11     /
          //          /            /            /
          //       q_4----------- x ---------- q_5
          //        /            /            /        z
          //       /    s_8     /    s_9     /         ^ y
          //      /            /            /          |/
          //    e_8----------q_2----------e_9          ---->x
          //
          //

          // loop over all coarse mesh cubes
          for(Index i(0); i < num_cubes; ++i)
          {
            // fetch coarse mesh quad-at-cube and edges-at-cube index vectors
            const IndexTupleTypeCQ& c_q = index_set_c_q[i];
            const IndexTupleTypeCE& c_e = index_set_c_e[i];

            // fetch fine mesh vertices-at-quad index vectors
            IndexTupleType& s_0 = index_set_out[offset + 12*i + 0];
            IndexTupleType& s_1 = index_set_out[offset + 12*i + 1];
            IndexTupleType& s_2 = index_set_out[offset + 12*i + 2];
            IndexTupleType& s_3 = index_set_out[offset + 12*i + 3];
            IndexTupleType& s_4 = index_set_out[offset + 12*i + 4];
            IndexTupleType& s_5 = index_set_out[offset + 12*i + 5];
            IndexTupleType& s_6 = index_set_out[offset + 12*i + 6];
            IndexTupleType& s_7 = index_set_out[offset + 12*i + 7];
            IndexTupleType& s_8 = index_set_out[offset + 12*i + 8];
            IndexTupleType& s_9 = index_set_out[offset + 12*i + 9];
            IndexTupleType& s_10 = index_set_out[offset + 12*i + 10];
            IndexTupleType& s_11 = index_set_out[offset + 12*i + 11];

            // calculate fine quad-vertex indices

            // x-normal
            s_0[0] = ioe + c_e[0]; // e_0
            s_0[1] = ioq + c_q[0]; // q_0
            s_0[2] = ioq + c_q[2]; // q_2
            s_0[3] = ioc + i;      // x

            s_1[0] = ioq + c_q[0]; // q_0
            s_1[1] = ioe + c_e[1]; // e_1
            s_1[2] = ioc + i;      // x
            s_1[3] = ioq + c_q[3]; // q_3

            s_2[0] = ioq + c_q[2]; // q_2
            s_2[1] = ioc + i;      // x
            s_2[2] = ioe + c_e[2]; // e_2
            s_2[3] = ioq + c_q[1]; // q_1

            s_3[0] = ioc + i;      // x
            s_3[1] = ioq + c_q[3]; // q_3
            s_3[2] = ioq + c_q[1]; // q_1
            s_3[3] = ioe + c_e[3]; // e_3

            // y-normal
            s_4[0] = ioe + c_e[4]; // e_4
            s_4[1] = ioq + c_q[0]; // q_0
            s_4[2] = ioq + c_q[4]; // q_4
            s_4[3] = ioc + i;      // x

            s_5[0] = ioq + c_q[0]; // q_0
            s_5[1] = ioe + c_e[5]; // e_5
            s_5[2] = ioc + i;      // x
            s_5[3] = ioq + c_q[5]; // q_5

            s_6[0] = ioq + c_q[4]; // q_4
            s_6[1] = ioc + i;      // x
            s_6[2] = ioe + c_e[6]; // e_6
            s_6[3] = ioq + c_q[1]; // q_1

            s_7[0] = ioc + i;      // x
            s_7[1] = ioq + c_q[5]; // q_5
            s_7[2] = ioq + c_q[1]; // q_1
            s_7[3] = ioe + c_e[7]; // e_7

            // z-normal
            s_8[0] = ioe + c_e[8]; // e_8
            s_8[1] = ioq + c_q[2]; // q_2
            s_8[2] = ioq + c_q[4]; // q_4
            s_8[3] = ioc + i;      // x

            s_9[0] = ioq + c_q[2]; // q_2
            s_9[1] = ioe + c_e[9]; // e_9
            s_9[2] = ioc + i;      // x
            s_9[3] = ioq + c_q[5]; // q_5

            s_10[0] = ioq + c_q[4]; // q_4
            s_10[1] = ioc + i;      // x
            s_10[2] = ioe + c_e[10]; // e_10
            s_10[3] = ioq + c_q[3]; // q_3

            s_11[0] = ioc + i;      // x
            s_11[1] = ioq + c_q[5]; // q_5
            s_11[2] = ioq + c_q[3]; // q_3
            s_11[3] = ioe + c_e[11]; // e_11
          }
          // return fine quad count
          return 12*num_cubes;
        }
      }; // StandardIndexRefiner<Hypercube<3>,2,0>

      /**
       * \brief StandardIndexRefiner implementation for Hypercube<3>: Vertices-at-Hexas
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Hypercube<3>, 3, 0>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 3; // second template parameter
        static constexpr int face_dim = 0; // third template parameter

        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioe = index_offsets[1];
          const Index ioq = index_offsets[2];
          const Index ioc = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types
          typedef IndexSet<6> IndexSetTypeCQ;
          typedef IndexSet<12> IndexSetTypeCE;
          typedef IndexSet<8> IndexSetTypeCV;

          // typedef index vector type
          typedef IndexSetTypeCQ::IndexTupleType IndexTupleTypeCQ;
          typedef IndexSetTypeCE::IndexTupleType IndexTupleTypeCE;
          typedef IndexSetTypeCV::IndexTupleType IndexTupleTypeCV;

          // fetch the index sets
          const IndexSetTypeCQ& index_set_c_q = index_set_holder_in.get_index_set<3,2>();
          const IndexSetTypeCE& index_set_c_e = index_set_holder_in.get_index_set<3,1>();
          const IndexSetTypeCV& index_set_c_v = index_set_holder_in.get_index_set<3,0>();

          // fetch number of coarse mesh cubes
          const Index num_cubes = index_set_c_q.get_num_entities();

          //
          //            v_6-----------e_3-----------v_7
          //             /|            /|            /|
          //            / |           / |           / |
          //           /  |          /  |          /  |
          //        e_6------------q_1-----------e_7  |
          //         /|   |        /|   |        /|   |
          //        / |  e_10-----/-|--q_3------/-|-e_11
          //       /  |  /|      /  |  /|      /  |  /|
          //     v_4-----------e_2-----------v_5  | / |
          //      |   |/  |     |   |/  |     |   |/  |
          //      |  q_4--------|---x---------|--q_5  |
          //      |  /|   |     |  /|   |     |  /|   |
          //      | / |  v_2----|-/-|--e_1----|-/-|--v_3
          //      |/  |  /      |/  |  /      |/  |  /
          //     e_8-----------q_2-----------e_9  | /
          //      |   |/        |   |/        |   |/
          //      |  e_4--------|--q_0--------|--e_5
          //      |  /          |  /          |  /
          //      | /           | /           | /
          //      |/            |/            |/
          //     v_0-----------e_0-----------v_1
          //

          // loop over all coarse mesh cubes
          for(Index i(0); i < num_cubes; ++i)
          {
            // fetch coarse mesh index vectors
            const IndexTupleTypeCQ& c_q = index_set_c_q[i];
            const IndexTupleTypeCE& c_e = index_set_c_e[i];
            const IndexTupleTypeCV& c_v = index_set_c_v[i];

            // fetch fine mesh vertices-at-cube index vectors
            IndexTupleType& c_0 = index_set_out[offset + 8*i + 0];
            IndexTupleType& c_1 = index_set_out[offset + 8*i + 1];
            IndexTupleType& c_2 = index_set_out[offset + 8*i + 2];
            IndexTupleType& c_3 = index_set_out[offset + 8*i + 3];
            IndexTupleType& c_4 = index_set_out[offset + 8*i + 4];
            IndexTupleType& c_5 = index_set_out[offset + 8*i + 5];
            IndexTupleType& c_6 = index_set_out[offset + 8*i + 6];
            IndexTupleType& c_7 = index_set_out[offset + 8*i + 7];

            // calculate fine cube-vertex indices

            c_0[0] = iov + c_v[0]; // v_0
            c_0[1] = ioe + c_e[0]; // e_0
            c_0[2] = ioe + c_e[4]; // e_4
            c_0[3] = ioq + c_q[0]; // q_0
            c_0[4] = ioe + c_e[8]; // e_8
            c_0[5] = ioq + c_q[2]; // q_2
            c_0[6] = ioq + c_q[4]; // q_4
            c_0[7] = ioc + i;      // x

            c_1[0] = ioe + c_e[0]; // e_0
            c_1[1] = iov + c_v[1]; // v_1
            c_1[2] = ioq + c_q[0]; // q_0
            c_1[3] = ioe + c_e[5]; // e_5
            c_1[4] = ioq + c_q[2]; // q_2
            c_1[5] = ioe + c_e[9]; // e_9
            c_1[6] = ioc + i;      // x
            c_1[7] = ioq + c_q[5]; // q_5

            c_2[0] = ioe + c_e[4]; // e_4
            c_2[1] = ioq + c_q[0]; // q_0
            c_2[2] = iov + c_v[2]; // v_2
            c_2[3] = ioe + c_e[1]; // e_1
            c_2[4] = ioq + c_q[4]; // q_4
            c_2[5] = ioc + i;      // x
            c_2[6] = ioe + c_e[10];// e_10
            c_2[7] = ioq + c_q[3]; // q_3

            c_3[0] = ioq + c_q[0]; // q_0
            c_3[1] = ioe + c_e[5]; // e_5
            c_3[2] = ioe + c_e[1]; // e_1
            c_3[3] = iov + c_v[3]; // v_3
            c_3[4] = ioc + i;      // x
            c_3[5] = ioq + c_q[5]; // q_5
            c_3[6] = ioq + c_q[3]; // q_3
            c_3[7] = ioe + c_e[11];// e_11

            c_4[0] = ioe + c_e[8]; // e_8
            c_4[1] = ioq + c_q[2]; // q_2
            c_4[2] = ioq + c_q[4]; // q_4
            c_4[3] = ioc + i;      // x
            c_4[4] = iov + c_v[4]; // v_4
            c_4[5] = ioe + c_e[2]; // e_2
            c_4[6] = ioe + c_e[6]; // e_6
            c_4[7] = ioq + c_q[1]; // q_1

            c_5[0] = ioq + c_q[2]; // q_2
            c_5[1] = ioe + c_e[9]; // e_9
            c_5[2] = ioc + i;      // x
            c_5[3] = ioq + c_q[5]; // q_5
            c_5[4] = ioe + c_e[2]; // e_2
            c_5[5] = iov + c_v[5]; // v_5
            c_5[6] = ioq + c_q[1]; // q_1
            c_5[7] = ioe + c_e[7]; // e_7

            c_6[0] = ioq + c_q[4]; // q_4
            c_6[1] = ioc + i;      // x
            c_6[2] = ioe + c_e[10];// e_10
            c_6[3] = ioq + c_q[3]; // q_3
            c_6[4] = ioe + c_e[6]; // e_6
            c_6[5] = ioq + c_q[1]; // q_1
            c_6[6] = iov + c_v[6]; // v_6
            c_6[7] = ioe + c_e[3]; // e_3

            c_7[0] = ioc + i;      // x
            c_7[1] = ioq + c_q[5]; // q_5
            c_7[2] = ioq + c_q[3]; // q_3
            c_7[3] = ioe + c_e[11];// e_11
            c_7[4] = ioq + c_q[1]; // q_1
            c_7[5] = ioe + c_e[7]; // e_7
            c_7[6] = ioe + c_e[3]; // e_3
            c_7[7] = iov + c_v[7]; // v_7
          }
          // return fine cube count
          return 8*num_cubes;
        }
      }; // StandardIndexRefiner<Hypercube<3>,3,0>

      /**
       * \brief StandardIndexRefiner implementation for Hypercube<3>: Edges-at-Quads
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Hypercube<3>, 2, 1>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 2; // second template parameter
        static constexpr int face_dim = 1; // third template parameter

        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch index offsets
          const Index ioq = index_offsets[2];
          const Index ioc = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types;
          typedef IndexSet<4> IndexSetTypeQV;
          typedef IndexSet<8> IndexSetTypeCV;
          typedef IndexSet<6> IndexSetTypeCQ;

          // typedef index vector type
          typedef IndexSetTypeCV::IndexTupleType IndexTupleTypeCV;
          typedef IndexSetTypeCQ::IndexTupleType IndexTupleTypeCQ;

          // fetch the index sets
          const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
          const IndexSetTypeCV& index_set_c_v = index_set_holder_in.get_index_set<3,0>();
          const IndexSetTypeCQ& index_set_c_q = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh cubes
          const Index num_cubes = index_set_c_q.get_num_entities();

          // typedef the sub-index mapping
          typedef Intern::SubIndexMapping<ShapeType, cell_dim, face_dim> SubIndexMappingType;

          // x-normal surfaces
          //                         e_3
          //                         /|
          //                        / |
          //                       /  |
          //                    q_1 s |
          //                     /| 3 |
          //                    / |   q_3
          //                   /  |  /|
          //                e_2 s | / |
          //                  | 2 |/  |
          //                  |   x s |
          //                  |  /| 1 |
          //                  | / |   e_1
          //                  |/  |  /
          //                q_2 s | /
          //                  | 0 |/
          //                  |  q_0
          //                  |  /                     z
          //                  | /                      ^ y
          //                  |/                       |/
          //                 e_0                       ---->x


          //  y-normal surfaces
          //
          //       e_6-----------q_1-----------e_7
          //        |             |             |
          //        |             |             |
          //        |     s_6     |     s_7     |
          //        |             |             |
          //        |             |             |
          //       q_4------------x------------q_5
          //        |             |             |
          //        |             |             |
          //        |     s_4     |     s_5     |      z
          //        |             |             |      ^ y
          //        |             |             |      |/
          //       e_4-----------q_0-----------e_5     ---->x


          // z-normal surfaces
          //
          //           e_10----------q_3----------e_11
          //            /            /            /
          //           /   s_10     /   s_11     /
          //          /            /            /
          //       q_4----------- x ---------- q_5
          //        /            /            /        z
          //       /    s_8     /    s_9     /         ^ y
          //      /            /            /          |/
          //    e_8----------q_2----------e_9          ---->x
          //
          //

          // loop over all coarse mesh cubes
          for(Index i(0); i < num_cubes; ++i)
          {
            // fetch coarse mesh vertices-at-quad and edges-at-quad index vectors
            const IndexTupleTypeCV& c_v = index_set_c_v[i];
            const IndexTupleTypeCQ& c_q = index_set_c_q[i];

            // create a sub-index mapping object
            SubIndexMappingType sim(c_v, c_q, index_set_q_v);

            // fetch fine mesh edges-at-quad index vectors
            IndexTupleType& q_0 = index_set_out[offset + 12*i + 0];
            IndexTupleType& q_1 = index_set_out[offset + 12*i + 1];
            IndexTupleType& q_2 = index_set_out[offset + 12*i + 2];
            IndexTupleType& q_3 = index_set_out[offset + 12*i + 3];
            IndexTupleType& q_4 = index_set_out[offset + 12*i + 4];
            IndexTupleType& q_5 = index_set_out[offset + 12*i + 5];
            IndexTupleType& q_6 = index_set_out[offset + 12*i + 6];
            IndexTupleType& q_7 = index_set_out[offset + 12*i + 7];
            IndexTupleType& q_8 = index_set_out[offset + 12*i + 8];
            IndexTupleType& q_9 = index_set_out[offset + 12*i + 9];
            IndexTupleType& q_10 = index_set_out[offset + 12*i + 10];
            IndexTupleType& q_11 = index_set_out[offset + 12*i + 11];

            // calculate fine edges-at-quad indices

            // x-normal
            q_0[0] = ioq + 4*c_q[0] + sim.map(0, 0);  // e_0-q_0
            q_0[1] = ioc + 6*i + 2;                   // q_2-x
            q_0[2] = ioq + 4*c_q[2] + sim.map(2, 0);  // e_0-q_2
            q_0[3] = ioc + 6*i + 0;                   // q_0-x

            q_1[0] = ioq + 4*c_q[0] + sim.map(0, 1);  // q_0-e_1
            q_1[1] = ioc + 6*i + 3;                   // x-q_3
            q_1[2] = ioc + 6*i + 0;                   // q_0-x
            q_1[3] = ioq + 4*c_q[3] + sim.map(3, 0);  // e_1-q_3

            q_2[0] = ioc + 6*i + 2;                   // q_2-x
            q_2[1] = ioq + 4*c_q[1] + sim.map(1, 0);  // e_2-q_1
            q_2[2] = ioq + 4*c_q[2] + sim.map(2, 1);  // q_2-e_2
            q_2[3] = ioc + 6*i + 1;                   // x-q_1

            q_3[0] = ioc + 6*i + 3;                   // x-q_3
            q_3[1] = ioq + 4*c_q[1] + sim.map(1, 1);  // q_1-e_3
            q_3[2] = ioc + 6*i + 1;                   // x-q_1
            q_3[3] = ioq + 4*c_q[3] + sim.map(3, 1);  // q_3-e_3

            // y-normal
            q_4[0] = ioq + 4*c_q[0] + sim.map(0, 2);  // e_4-q_0
            q_4[1] = ioc + 6*i + 4;                   // q_4-x
            q_4[2] = ioq + 4*c_q[4] + sim.map(4, 0);  // e_4-q_4
            q_4[3] = ioc + 6*i + 0;                   // q_0-x

            q_5[0] = ioq + 4*c_q[0] + sim.map(0, 3);  // q_0-e_5
            q_5[1] = ioc + 6*i + 5;                   // x-q_5
            q_5[2] = ioc + 6*i + 0;                   // q_0-x
            q_5[3] = ioq + 4*c_q[5] + sim.map(5, 0);  // e_5-q_5

            q_6[0] = ioc + 6*i + 4;                   // q_4-x
            q_6[1] = ioq + 4*c_q[1] + sim.map(1, 2);  // e_6-q_1
            q_6[2] = ioq + 4*c_q[4] + sim.map(4, 1);  // q_4-e_6
            q_6[3] = ioc + 6*i + 1;                   // x-q_1

            q_7[0] = ioc + 6*i + 5;                   // x-q_5
            q_7[1] = ioq + 4*c_q[1] + sim.map(1, 3);  // q_1-e_7
            q_7[2] = ioc + 6*i + 1;                   // x-q_1
            q_7[3] = ioq + 4*c_q[5] + sim.map(5, 1);  // q_5-e_7

            // z-normal
            q_8[0] = ioq + 4*c_q[2] + sim.map(2, 2);  // e_8-q_2
            q_8[1] = ioc + 6*i + 4;                   // q_4-x
            q_8[2] = ioq + 4*c_q[4] + sim.map(4, 2);  // e_8-q_4
            q_8[3] = ioc + 6*i + 2;                   // q_2-x

            q_9[0] = ioq + 4*c_q[2] + sim.map(2, 3);  // q_2-e_9
            q_9[1] = ioc + 6*i + 5;                   // x-q_5
            q_9[2] = ioc + 6*i + 2;                   // q_2-x
            q_9[3] = ioq + 4*c_q[5] + sim.map(5, 2);  // e_9-q_5

            q_10[0] = ioc + 6*i + 4;                   // q_4-x
            q_10[1] = ioq + 4*c_q[3] + sim.map(3, 2);  // e_10-q_3
            q_10[2] = ioq + 4*c_q[4] + sim.map(4, 3);  // q_4-e_10
            q_10[3] = ioc + 6*i + 3;                   // x-q_3

            q_11[0] = ioc + 6*i + 5;                   // x-q_5
            q_11[1] = ioq + 4*c_q[3] + sim.map(3, 3);  // q_3-e_11
            q_11[2] = ioc + 6*i + 3;                   // x-q_3
            q_11[3] = ioq + 4*c_q[5] + sim.map(5, 3);  // q_5-e_11
          }
          // return fine quad count
          return 12*num_cubes;
        }
      }; // StandardIndexRefiner<Hypercube<3>,2,1>

      /**
       * \brief StandardIndexRefiner implementation for Hypercube<3>: Edges-at-Hexas
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Hypercube<3>, 3, 1>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 3; // second template parameter
        static constexpr int face_dim = 1; // third template parameter

        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch edge index offsets
          const Index ioe = index_offsets[1];
          const Index ioq = index_offsets[2];
          const Index ioc = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types
          typedef IndexSet<2> IndexSetTypeEV;
          typedef IndexSet<4> IndexSetTypeQV;
          typedef IndexSet<8> IndexSetTypeCV;
          typedef IndexSet<12> IndexSetTypeCE;
          typedef IndexSet<6> IndexSetTypeCQ;

          // typedef index vector type
          typedef IndexSetTypeCV::IndexTupleType IndexTupleTypeCV;
          typedef IndexSetTypeCE::IndexTupleType IndexTupleTypeCE;
          typedef IndexSetTypeCQ::IndexTupleType IndexTupleTypeCQ;

          // fetch the index sets
          const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();
          const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
          const IndexSetTypeCV& index_set_c_v = index_set_holder_in.get_index_set<3,0>();
          const IndexSetTypeCE& index_set_c_e = index_set_holder_in.get_index_set<3,1>();
          const IndexSetTypeCQ& index_set_c_q = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh cubes
          const Index num_cubes = index_set_c_q.get_num_entities();

          // typedef the sub-index mapping
          typedef Intern::SubIndexMapping<ShapeType, 2, 1> SubIndexMappingType;
          typedef Intern::SubIndexMapping<ShapeType, face_dim, 0> EdgeIndexMappingType;

          //
          //            v_6-----------e_3-----------v_7
          //             /|            /|            /|
          //            / |           / |           / |
          //           /  |          /  |          /  |
          //        e_6------------q_1-----------e_7  |
          //         /|   |        /|   |        /|   |
          //        / |  e_10-----/-|--q_3------/-|-e_11
          //       /  |  /|      /  |  /|      /  |  /|
          //     v_4-----------e_2-----------v_5  | / |
          //      |   |/  |     |   |/  |     |   |/  |
          //      |  q_4--------|---x---------|--q_5  |
          //      |  /|   |     |  /|   |     |  /|   |
          //      | / |  v_2----|-/-|--e_1----|-/-|--v_3
          //      |/  |  /      |/  |  /      |/  |  /
          //     e_8-----------q_2-----------e_9  | /
          //      |   |/        |   |/        |   |/
          //      |  e_4--------|--q_0--------|--e_5
          //      |  /          |  /          |  /
          //      | /           | /           | /
          //      |/            |/            |/
          //     v_0-----------e_0-----------v_1
          //

          // loop over all coarse mesh cubes
          for(Index i(0); i < num_cubes; ++i)
          {
            // fetch coarse mesh index vectors
            const IndexTupleTypeCV& c_v = index_set_c_v[i];
            const IndexTupleTypeCE& c_e = index_set_c_e[i];
            const IndexTupleTypeCQ& c_q = index_set_c_q[i];

            // create sub-index mapping objects
            SubIndexMappingType sim(c_v, c_q, index_set_q_v);
            EdgeIndexMappingType edgesim(c_v, c_e, index_set_e_v);

            // fetch fine mesh edges-at-cube index vectors
            IndexTupleType& c_0 = index_set_out[offset + 8*i + 0];
            IndexTupleType& c_1 = index_set_out[offset + 8*i + 1];
            IndexTupleType& c_2 = index_set_out[offset + 8*i + 2];
            IndexTupleType& c_3 = index_set_out[offset + 8*i + 3];
            IndexTupleType& c_4 = index_set_out[offset + 8*i + 4];
            IndexTupleType& c_5 = index_set_out[offset + 8*i + 5];
            IndexTupleType& c_6 = index_set_out[offset + 8*i + 6];
            IndexTupleType& c_7 = index_set_out[offset + 8*i + 7];

            // calculate fine edges-at-cube indices

            // v_0-cube
            c_0[0] = ioe + 2*c_e[0] + edgesim.map(0, 0);    // v_0-e_0
            c_0[1] = ioq + 4*c_q[0] + sim.map(0, 2);        // e_4-q_0
            c_0[2] = ioq + 4*c_q[2] + sim.map(2, 2);        // e_8-q_2
            c_0[3] = ioc + 6*i + 4;                         // q_4-x
            c_0[4] = ioe + 2*c_e[4] + edgesim.map(4, 0);    // v_0-e_4
            c_0[5] = ioq + 4*c_q[0] + sim.map(0, 0);        // e_0-q_0
            c_0[6] = ioq + 4*c_q[4] + sim.map(4, 2);        // e_8-q_4
            c_0[7] = ioc + 6*i + 2;                         // q_2-x
            c_0[8] = ioe + 2*c_e[8] + edgesim.map(8, 0);    // v_0-e_8
            c_0[9] = ioq + 4*c_q[2] + sim.map(2, 0);        // e_0-q_2
            c_0[10] = ioq + 4*c_q[4] + sim.map(4, 0);       // e_4-q_4
            c_0[11] = ioc + 6*i + 0;                        // q_0-x

            // v_1-cube
            c_1[0] = ioe + 2*c_e[0] + edgesim.map(0,1);     // e_0-v_1
            c_1[1] = ioq + 4*c_q[0] + sim.map(0, 3);        // q_0-e_5
            c_1[2] = ioq + 4*c_q[2] + sim.map(2,3);         // q_2-e_9
            c_1[3] = ioc + 6*i + 5;                         // x-q_5
            c_1[4] = ioq + 4*c_q[0] + sim.map(0, 0);        // e_0-q_0
            c_1[5] = ioe + 2*c_e[5] + edgesim.map(5, 0);    // v_1-e_5
            c_1[6] = ioc + 6*i + 2;                         // q_2-x
            c_1[7] = ioq + 4*c_q[5] + sim.map(5, 2);        // e_9-q_5
            c_1[8] = ioq + 4*c_q[2] + sim.map(2,0);         // e_0-q_2
            c_1[9] = ioe + 2*c_e[9] + edgesim.map(9, 0);    // v_1-e_9
            c_1[10] = ioc + 6*i + 0;                        // q_0-x
            c_1[11] = ioq + 4*c_q[5] + sim.map(5, 0);       // e_5-q_5

            // v_2-cube
            c_2[0] = ioq + 4*c_q[0] + sim.map(0, 2);        // e_4-q_0
            c_2[1] = ioe + 2*c_e[1] + edgesim.map(1, 0);    // v_2-e_1
            c_2[2] = ioc + 6*i + 4;                         // q_4-x
            c_2[3] = ioq + 4*c_q[3] + sim.map(3, 2);        // e_10-q_3
            c_2[4] = ioe + 2*c_e[4] + edgesim.map(4, 1);    // e_4-v_2
            c_2[5] = ioq + 4*c_q[0] + sim.map(0, 1);        // q_0-e_1
            c_2[6] = ioq + 4*c_q[4] + sim.map(4, 3);        // q_4-e_10
            c_2[7] = ioc + 6*i + 3;                         // x-q_3
            c_2[8] = ioq + 4*c_q[4] + sim.map(4, 0);        // e_4-q_4
            c_2[9] = ioc + 6*i + 0;                         // q_0-x
            c_2[10] = ioe + 2*c_e[10] + edgesim.map(10, 0); // v_2-e_10
            c_2[11] = ioq + 4*c_q[3] + sim.map(3, 0);       // e_1-q_3

            // v_3-cube
            c_3[0] = ioq + 4*c_q[0] + sim.map(0, 3);        // q_0-e_5
            c_3[1] = ioe + 2*c_e[1] + edgesim.map(1, 1);    // e_1-v_3
            c_3[2] = ioc + 6*i + 5;                         // x-q_5
            c_3[3] = ioq + 4*c_q[3] + sim.map(3, 3);        // q_3-e_11
            c_3[4] = ioq + 4*c_q[0] + sim.map(0, 1);        // q_0-e_1
            c_3[5] = ioe + 2*c_e[5] + edgesim.map(5, 1);    // e_5-v_3
            c_3[6] = ioc + 6*i + 3;                         // x-q_3
            c_3[7] = ioq + 4*c_q[5] + sim.map(5, 3);        // q_5-e_11
            c_3[8] = ioc + 6*i + 0;                         // q_0-x
            c_3[9] = ioq + 4*c_q[5] + sim.map(5, 0);        // e_5-q_5
            c_3[10] = ioq + 4*c_q[3] + sim.map(3, 0);       // e_1-q_3
            c_3[11] = ioe + 2*c_e[11] + edgesim.map(11, 0); // v_3-e_11

            // v_4-cube
            c_4[0] = ioq + 4*c_q[2] + sim.map(2, 2);        // e_8-q_2
            c_4[1] = ioc + 6*i + 4;                         // q_4-x
            c_4[2] = ioe + 2*c_e[2] + edgesim.map(2, 0);    // v_4-e_2
            c_4[3] = ioq + 4*c_q[1] + sim.map(1, 2);        // e_6-q_1
            c_4[4] = ioq + 4*c_q[4] + sim.map(4, 2);        // e_8-q_4
            c_4[5] = ioc + 6*i + 2;                         // q_2-x
            c_4[6] = ioe + 2*c_e[6] + edgesim.map(6, 0);    // v_4-e_6
            c_4[7] = ioq + 4*c_q[1] + sim.map(1, 0);        // e_2-q_1
            c_4[8] = ioe + 2*c_e[8] + edgesim.map(8, 1);    // e_8-v_4
            c_4[9] = ioq + 4*c_q[2] + sim.map(2, 1);        // q_2-e_2
            c_4[10] = ioq + 4*c_q[4] + sim.map(4, 1);       // q_4-e_6
            c_4[11] = ioc + 6*i + 1;                         // x-q_1

            // v_5-cube
            c_5[0] = ioq + 4*c_q[2] + sim.map(2, 3);        // q_2-e_9
            c_5[1] = ioc + 6*i + 5;                         // x-q_5
            c_5[2] = ioe + 2*c_e[2] + edgesim.map(2, 1);    // e_2-v-5
            c_5[3] = ioq + 4*c_q[1] + sim.map(1, 3);        // q_1-e_7
            c_5[4] = ioc + 6*i + 2;                         // q_2-x
            c_5[5] = ioq + 4*c_q[5] + sim.map(5, 2);        // e_9-q_5
            c_5[6] = ioq + 4*c_q[1] + sim.map(1, 0);        // e_2-q_1
            c_5[7] = ioe + 2*c_e[7] + edgesim.map(7, 0);    // v_5-e_7
            c_5[8] = ioq + 4*c_q[2] + sim.map(2, 1);        // q_2-e_2
            c_5[9] = ioe + 2*c_e[9] + edgesim.map(9, 1);    // e_9-v_5
            c_5[10] = ioc + 6*i + 1;                        // x-q_1
            c_5[11] = ioq + 4*c_q[5] + sim.map(5, 1);       // q_5-e_7

            // v_6-cube
            c_6[0] = ioc + 6*i + 4;                         // q_4-x
            c_6[1] = ioq + 4*c_q[3] + sim.map(3, 2);        // e_10-q_3
            c_6[2] = ioq + 4*c_q[1] + sim.map(1, 2);        // e_6-q_1
            c_6[3] = ioe + 2*c_e[3] + edgesim.map(3, 0);    // v_6-e_3
            c_6[4] = ioq + 4*c_q[4] + sim.map(4, 3);        // q_4-e_10
            c_6[5] = ioc + 6*i + 3;                         // x-q_3
            c_6[6] = ioe + 2*c_e[6] + edgesim.map(6, 1);    // e_6-v_6
            c_6[7] = ioq + 4*c_q[1] + sim.map(1, 1);        // q_1-e_3
            c_6[8] = ioq + 4*c_q[4] + sim.map(4, 1);        // q_4-e_6
            c_6[9] = ioc + 6*i + 1;                         // x-q_1
            c_6[10] = ioe + 2*c_e[10] + edgesim.map(10, 1); // e_10-v_6
            c_6[11] = ioq + 4*c_q[3] + sim.map(3, 1);       // q_3-e_3

            // v_7-cube
            c_7[0] = ioc + 6*i + 5;                         // x-q_5
            c_7[1] = ioq + 4*c_q[3] + sim.map(3, 3);        // q_3-e_11
            c_7[2] = ioq + 4*c_q[1] + sim.map(1, 3);        // q_1-e_7
            c_7[3] = ioe + 2*c_e[3] + edgesim.map(3, 1);    // e_3-v_7
            c_7[4] = ioc + 6*i + 3;                         // x-q_3
            c_7[5] = ioq + 4*c_q[5] + sim.map(5, 3);        // q_5-e_11
            c_7[6] = ioq + 4*c_q[1] + sim.map(1, 1);        // q_1-e_3
            c_7[7] = ioe + 2*c_e[7] + edgesim.map(7, 1);    // e_7-v_7
            c_7[8] = ioc + 6*i + 1;                         // x-q_1
            c_7[9] = ioq + 4*c_q[5] + sim.map(5, 1);        // q_5-e_7
            c_7[10] = ioq + 4*c_q[3] + sim.map(3, 1);       // q_3-e_3
            c_7[11] = ioe + 2*c_e[11] + edgesim.map(11, 1); // e_11-v_7
          }
          // return fine cube count
          return 8*num_cubes;
        }
      }; // StandardIndexRefiner<Hypercube<3>,3,1>

      /**
       * \brief StandardIndexRefiner implementation for Hypercube<3>: Quads-at-Hexas
       *
       * \author Constantin Christof
       */
      template<>
      struct StandardIndexRefiner<Shape::Hypercube<3>, 3, 2>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 3; // second template parameter
        static constexpr int face_dim = 2; // third template parameter

        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch quad index offsets
          const Index ioq = index_offsets[2];
          const Index ioc = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexTupleType IndexTupleType;

          // typedef index set types
          typedef IndexSet<4> IndexSetTypeQV;
          typedef IndexSet<8> IndexSetTypeCV;
          typedef IndexSet<6> IndexSetTypeCQ;

          // typedef index vector type
          typedef IndexSetTypeCV::IndexTupleType IndexTupleTypeCV;
          typedef IndexSetTypeCQ::IndexTupleType IndexTupleTypeCQ;

          // fetch the index sets
          const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
          const IndexSetTypeCV& index_set_c_v = index_set_holder_in.get_index_set<3,0>();
          const IndexSetTypeCQ& index_set_c_q = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh cubes
          const Index num_cubes = index_set_c_q.get_num_entities();

          // typedef the sub-index mapping
          typedef Intern::SubIndexMapping<ShapeType, face_dim, 0> SubIndexMappingType;

          //
          //            v_6-----------e_3-----------v_7
          //             /|            /|            /|
          //            / |           / |           / |
          //           /  |          /  |          /  |
          //        e_6------------q_1-----------e_7  |
          //         /|   |        /|   |        /|   |
          //        / |  e_10-----/-|--q_3------/-|-e_11
          //       /  |  /|      /  |  /|      /  |  /|
          //     v_4-----------e_2-----------v_5  | / |
          //      |   |/  |     |   |/  |     |   |/  |
          //      |  q_4--------|---x---------|--q_5  |
          //      |  /|   |     |  /|   |     |  /|   |
          //      | / |  v_2----|-/-|--e_1----|-/-|--v_3
          //      |/  |  /      |/  |  /      |/  |  /
          //     e_8-----------q_2-----------e_9  | /
          //      |   |/        |   |/        |   |/
          //      |  e_4--------|--q_0--------|--e_5
          //      |  /          |  /          |  /
          //      | /           | /           | /
          //      |/            |/            |/
          //     v_0-----------e_0-----------v_1
          //

          // loop over all coarse mesh cubes
          for(Index i(0); i < num_cubes; ++i)
          {
            // fetch coarse mesh vertices-at-cube and quads-at-cube index vectors
            const IndexTupleTypeCV& c_v = index_set_c_v[i];
            const IndexTupleTypeCQ& c_q = index_set_c_q[i];

            // create a sub-index mapping object
            SubIndexMappingType sim(c_v, c_q, index_set_q_v);

            // fetch fine mesh quads-at-cube index vectors
            IndexTupleType& c_0 = index_set_out[offset + 8*i + 0];
            IndexTupleType& c_1 = index_set_out[offset + 8*i + 1];
            IndexTupleType& c_2 = index_set_out[offset + 8*i + 2];
            IndexTupleType& c_3 = index_set_out[offset + 8*i + 3];
            IndexTupleType& c_4 = index_set_out[offset + 8*i + 4];
            IndexTupleType& c_5 = index_set_out[offset + 8*i + 5];
            IndexTupleType& c_6 = index_set_out[offset + 8*i + 6];
            IndexTupleType& c_7 = index_set_out[offset + 8*i + 7];

            // calculate fine quad-at-cube indices

            // v_0-cube
            c_0[0] = ioq + 4*c_q[0] + sim.map(0, 0);  // v_0-e_0-e_4-q_0
            c_0[1] = ioc + 12*i + 8;                  // e_8-q_2-q_4-x
            c_0[2] = ioq + 4*c_q[2] + sim.map(2, 0);  // v_0-e_0-e_8-q_2
            c_0[3] = ioc + 12*i + 4;                  // e_4-q_0-q_4-x
            c_0[4] = ioq + 4*c_q[4] + sim.map(4, 0);  // v_0-e_4-e_8-q_4
            c_0[5] = ioc + 12*i + 0;                  // e_0-q_0-q_2-x


            // v_1-cube
            c_1[0] = ioq + 4*c_q[0] + sim.map(0, 1);  // e_0-v_1-q_0-e_5
            c_1[1] = ioc + 12*i + 9;                  // q_2-e_9-x-q_5
            c_1[2] = ioq + 4*c_q[2] + sim.map(2, 1);  // e_0-v_1-q_2-e_9
            c_1[3] = ioc + 12*i + 5;                  // q_0-e_5-x-q_5
            c_1[4] = ioc + 12*i + 0;                  // e_0-q_0-q_2-x
            c_1[5] = ioq + 4*c_q[5] + sim.map(5, 0);  // v_1-e_5-e_9-q_5

            // v_2-cube
            c_2[0] = ioq + 4*c_q[0] + sim.map(0, 2);  // e_4-q_0-v_2-e_1
            c_2[1] = ioc + 12*i + 10;                 // q_4-x-e_10-q_3
            c_2[2] = ioc + 12*i + 4;                  // e_4-q_0-q_4-x
            c_2[3] = ioq + 4*c_q[3] + sim.map(3, 0);  // v_2-e_1-e_10-q_3
            c_2[4] = ioq + 4*c_q[4] + sim.map(4, 1);  // e_4-v_2-q_4-e_10
            c_2[5] = ioc + 12*i + 1;                  // q_0-e_1-x-q_3

            // v_3-cube
            c_3[0] = ioq + 4*c_q[0] + sim.map(0, 3);  // q_0-e_5-e_1-v_3
            c_3[1] = ioc + 12*i + 11;                 // x-q_5-q_3-e_11
            c_3[2] = ioc + 12*i + 5;                  // q_0-e_5-x-q_5
            c_3[3] = ioq + 4*c_q[3] + sim.map(3, 1);  // e_1-v_3-q_3-e_11
            c_3[4] = ioc + 12*i + 1;                  // q_0-e_1-x-q_3
            c_3[5] = ioq + 4*c_q[5] + sim.map(5, 1);  // e_5-v_3-q_5-e_11

            // v_4-cube
            c_4[0] = ioc + 12*i + 8;                  // e_8-q_2-q_4-x
            c_4[1] = ioq + 4*c_q[1] + sim.map(1, 0);  // v_4-e_2-e_6-q_1
            c_4[2] = ioq + 4*c_q[2] + sim.map(2, 2);  // e_8-q_2-v_4-e_2
            c_4[3] = ioc + 12*i + 6;                  // q_4-x-e_6-q_1
            c_4[4] = ioq + 4*c_q[4] + sim.map(4, 2);  // e_8-q_4-v_4-e_6
            c_4[5] = ioc + 12*i + 2;                  // q_2-x-e_2-q_1

            // v_5-cube
            c_5[0] = ioc + 12*i + 9;                  // q_2-e_9-x-q_5
            c_5[1] = ioq + 4*c_q[1] + sim.map(1, 1);  // e_2-v_5-q_1-e_7
            c_5[2] = ioq + 4*c_q[2] + sim.map(2, 3);  // q_2-e_9-e_2-v_5
            c_5[3] = ioc + 12*i + 7;                  // x-q_5-q_1-e_7
            c_5[4] = ioc + 12*i + 2;                  // q_2-x-e_2-q_1
            c_5[5] = ioq + 4*c_q[5] + sim.map(5, 2);  // e_9-q_5-v_5-e_7

            // v_6-cube
            c_6[0] = ioc + 12*i + 10;                 // q_4-x-e_10-q_3
            c_6[1] = ioq + 4*c_q[1] + sim.map(1, 2);  // e_6-q_1-v_6-e_3
            c_6[2] = ioc + 12*i + 6;                  // q_4-x-e_6-q_1
            c_6[3] = ioq + 4*c_q[3] + sim.map(3, 2);  // e_10-q_3-v_6-e_3
            c_6[4] = ioq + 4*c_q[4] + sim.map(4, 3);  // q_4-e_10-e_6-v_6
            c_6[5] = ioc + 12*i + 3;                  // x-q_3-q_1-e_3

            // v_7-cube
            c_7[0] = ioc + 12*i + 11;                 // x-q_5-q_3-e_11
            c_7[1] = ioq + 4*c_q[1] + sim.map(1, 3);  // q_1-e_7-e_3-v_7
            c_7[2] = ioc + 12*i + 7;                  // x-q_5-q_1-e_7
            c_7[3] = ioq + 4*c_q[3] + sim.map(3, 3);  // q_3-e_11-e_3-v_7
            c_7[4] = ioc + 12*i + 3;                  // x-q_3-q_1-e_3
            c_7[5] = ioq + 4*c_q[5] + sim.map(5, 3);  // q_5-e_11-e_7-v_7
          }
          // return fine cube count
          return 8*num_cubes;
        }
      }; // StandardIndexRefiner<Hypercube<3>,3,2>

      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */

      template<
        typename Shape_,
        int cell_dim_,
        int face_dim_,
        // The following "dummy" argument is necessary for partial specialisation;
        // shape_dim_ *must* coincide with Shape_::dimension !!!
        int shape_dim_ = Shape_::dimension>
      struct IndexRefineShapeWrapper
      {
        static_assert(shape_dim_ == Shape_::dimension, "invalid shape dimension");
        // the case shape_dim_ = cell_dim_ is specialised below
        static_assert(shape_dim_ > cell_dim_, "invalid cell dimension");
        static_assert(cell_dim_ > face_dim_, "invalid face dimension");
        static_assert(face_dim_ >= 0, "invalid face_dimension");

        typedef Shape_ ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim_>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static String name()
        {
          return "IndexRefineShapeWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ","
            + stringify(face_dim_) + "," + stringify(shape_dim_) + ">";
        }

        static Index refine(
          IndexSetType& index_set_out,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          typedef typename Shape::FaceTraits<ShapeType, shape_dim_ - 1>::ShapeType FacetType;

          // refine facets
          Index offset = IndexRefineShapeWrapper<FacetType, cell_dim_, face_dim_>
            ::refine(index_set_out, index_offsets, index_set_holder_in);

          // call index refiner and return new offset
          Index num_faces = StandardIndexRefiner<Shape_, cell_dim_, face_dim_>
            ::refine(index_set_out, offset, index_offsets, index_set_holder_in);

          // validate number of created indices
          Index num_shapes = index_set_holder_in.template get_index_set<shape_dim_, 0>().get_num_entities();
          XASSERTM(num_faces == (StandardRefinementTraits<ShapeType, cell_dim_>::count * num_shapes),
            "IndexRefiner output does not match StdRefTraits prediction");

          // calculate new offset and return
          return offset + num_faces;
        }
      };

      template<
        typename Shape_,
        int cell_dim_,
        int face_dim_>
      struct IndexRefineShapeWrapper<Shape_, cell_dim_, face_dim_, cell_dim_>
      {
        static_assert(cell_dim_ == Shape_::dimension, "invalid shape dimension");
        static_assert(cell_dim_ > face_dim_, "invalid face dimension");
        static_assert(face_dim_ >= 0, "invalid face_dimension");

        typedef Shape_ ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim_>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static String name()
        {
          return "IndexRefineShapeWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ","
            + stringify(face_dim_) + "," + stringify(cell_dim_) + ">";
        }

        static Index refine(
          IndexSetType& index_set_out,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // call index refiner
          Index num_faces = StandardIndexRefiner<Shape_, cell_dim_, face_dim_>
            ::refine(index_set_out, 0, index_offsets, index_set_holder_in);

          // validate number of created indices
          Index num_shapes = index_set_holder_in.template get_index_set<cell_dim_, 0>().get_num_entities();
          XASSERTM(num_faces == (StandardRefinementTraits<ShapeType, cell_dim_>::count * num_shapes),
            "IndexRefiner output does not match StdRefTraits prediction");

          // return offset
          return num_faces;
        }
      };

      /* ********************************************************************************************************* */

      template<
        typename Shape_,
        int cell_dim_,
        int face_dim_ = cell_dim_ - 1>
      struct IndexRefineFaceWrapper
      {
        // the case face_dim_ = 0 is handled by the partial specialisation below
        static_assert(face_dim_ > 0, "invalid face dimension");
        static_assert(cell_dim_ > face_dim_, "invalid cell dimension");
        static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");

        typedef Shape_ ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, face_dim_>::count> IndexSetType;
        typedef IndexSetWrapper<CellType, face_dim_> IndexSetWrapperType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static String name()
        {
          return "IndexRefineFaceWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ","
            + stringify(face_dim_) + ">";
        }

        static void refine(
          IndexSetWrapperType& index_set_wrapper_out,
          const Index num_entities[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // recursive call of IndexRefineFaceWrapper
          IndexRefineFaceWrapper<Shape_, cell_dim_, face_dim_ - 1>
            ::refine(index_set_wrapper_out, num_entities, index_set_holder_in);

          // fetch output index set
          IndexSetType& index_set_out = index_set_wrapper_out.template get_index_set<face_dim_>();

          // calculate index offsets
          Index index_offsets[Shape_::dimension+1];
          EntityCounter<StandardRefinementTraits, Shape_, face_dim_>::offset(index_offsets, num_entities);

          // call refinement shape wrapper
          IndexRefineShapeWrapper<Shape_, cell_dim_, face_dim_>
            ::refine(index_set_out, index_offsets, index_set_holder_in);
        }
      };

      template<
        typename Shape_,
        int cell_dim_>
      struct IndexRefineFaceWrapper<Shape_, cell_dim_, 0>
      {
        static_assert(cell_dim_ > 0, "invalid cell dimension");
        static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");

        typedef Shape_ ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetWrapper<CellType, 0> IndexSetWrapperType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static String name()
        {
          return "IndexRefineFaceWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ",0>";
        }

        static void refine(
          IndexSetWrapperType& index_set_wrapper_out,
          const Index num_entities[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch output index set
          IndexSetType& index_set_out = index_set_wrapper_out.template get_index_set<0>();

          // calculate index offsets
          Index index_offsets[Shape_::dimension+1];
          EntityCounter<StandardRefinementTraits, Shape_, 0>::offset(index_offsets, num_entities);

          // call refinement shape wrapper
          IndexRefineShapeWrapper<Shape_, cell_dim_, 0>
            ::refine(index_set_out, index_offsets, index_set_holder_in);
        }
      };

      /* ********************************************************************************************************* */

      template<
        typename Shape_,
        int cell_dim_ = Shape_::dimension>
      struct IndexRefineWrapper
      {
        // the case cell_dim_ = 1 is handled by the partial specialisation below
        static_assert(cell_dim_ > 1, "invalid cell dimension");
        static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");

        typedef Shape_ ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef IndexSetWrapper<CellType> IndexSetWrapperType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static String name()
        {
          return "IndexRefineWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ">";
        }

        static void refine(
          IndexSetHolderType& index_set_holder_out,
          const Index num_entities[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // recursive call of IndexRefineWrapper
          IndexRefineWrapper<Shape_, cell_dim_ - 1>
            ::refine(index_set_holder_out, num_entities, index_set_holder_in);

          // fetch output index set wrapper
          IndexSetWrapperType& index_set_wrapper_out = index_set_holder_out
            .template get_index_set_wrapper<cell_dim_>();

          // call face wrapper
          IndexRefineFaceWrapper<Shape_, cell_dim_>
            ::refine(index_set_wrapper_out, num_entities, index_set_holder_in);
        }
      };

      template<typename Shape_>
      struct IndexRefineWrapper<Shape_, 1>
      {
        static_assert(1 <= Shape_::dimension, "invalid shape dimension");

        typedef Shape_ ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, 1>::ShapeType CellType;
        typedef IndexSetWrapper<CellType> IndexSetWrapperType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;

        static String name()
        {
          return "IndexRefineWrapper<" + Shape_::name() + ",1>";
        }

        static void refine(
          IndexSetHolderType& index_set_holder_out,
          const Index num_entities[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch output index set wrapper
          IndexSetWrapperType& index_set_wrapper_out = index_set_holder_out.template get_index_set_wrapper<1>();

          // call face wrapper
          IndexRefineFaceWrapper<Shape_, 1>::refine(index_set_wrapper_out, num_entities, index_set_holder_in);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT


#endif // KERNEL_GEOMETRY_INTERN_STANDARD_INDEX_REFINER_HPP
