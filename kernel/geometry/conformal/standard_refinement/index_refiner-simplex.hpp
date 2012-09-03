#pragma once
#ifndef KERNEL_GEOMETRY_CONF_STD_REF_IDX_REF_SIMPLEX_HPP
#define KERNEL_GEOMETRY_CONF_STD_REF_IDX_REF_SIMPLEX_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal/standard_refinement/index_refiner.hpp>
#include <kernel/geometry/congruency/sub_index_mapping.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Conformal
    {
      /// \cond internal
      namespace StandardRefinement
      {
        /**
         * \brief IndexRefiner implementation for Simplex<1>: Vertices-at-Edges
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Simplex<1>, 1, 0>
        {
          enum
          {
            shape_dim = 1,
            cell_dim = 1, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set type; an edge has 2 vertices, so we need an IndexSet<2>
            typedef IndexSet<2> IndexSetTypeEV;

            // typedef index vector reference type
            typedef IndexSetTypeEV::ConstIndexVectorReference ConstIndexVectorReferenceEV;

            // fetch the vertices-at-edge index set
            const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();

            // fetch number of coarse mesh edges
            const Index num_edges = index_set_e_v.get_num_entities();

            // loop over all coarse mesh edges
            for(Index i(0); i < num_edges; ++i)
            {
              // fetch coarse mesh vertices-at-edge index vector
              ConstIndexVectorReferenceEV e_v = index_set_e_v[i];

              // fetch fine mesh vertices-at-edge index vectors
              IndexVectorReference idx0 = index_set_out[offset + 2*i + 0];
              IndexVectorReference idx1 = index_set_out[offset + 2*i + 1];

              // calculate fine edge-vertex indices
              idx0[0] = iov + e_v[0];
              idx0[1] = ioe + i;
              idx1[0] = ioe + i;
              idx1[1] = iov + e_v[1];
            }

            // return fine edge count
            return 2*num_edges;
          }
        }; // IndexRefiner<Simplex<1>,1,0>

        /* ************************************************************************************* */
        /* ************************************************************************************* */
        /* ************************************************************************************* */

        /**
         * \brief IndexRefiner implementation for Simplex<2>: Vertices-at-Edges
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Simplex<2>, 1, 0>
        {
          enum
          {
            shape_dim = 2,
            cell_dim = 1, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set type; a triangle has 3 edges, so we need an IndexSet<3>
            typedef IndexSet<3> IndexSetTypeSE;

            // typedef index vector type
            typedef IndexSetTypeSE::ConstIndexVectorReference ConstIndexVectorReferenceSE;

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
              ConstIndexVectorReferenceSE s_e = index_set_s_e[i];

              // fetch fine mesh vertices-at-edge index vectors
              IndexVectorReference e_0 = index_set_out[offset + 3*i + 0];
              IndexVectorReference e_1 = index_set_out[offset + 3*i + 1];
              IndexVectorReference e_2 = index_set_out[offset + 3*i + 2];

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
        }; // IndexRefiner<Simplex<2>,1,0>

        /**
         * \brief IndexRefiner implementation for Simplex<2>: Vertices-at-Tris
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Simplex<2>, 2, 0>
        {
          enum
          {
            shape_dim = 2,
            cell_dim = 2, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types; a triangle has 3 vertices and 3 edges, so both need an IndexSet<3>
            typedef IndexSet<3> IndexSetTypeSV;
            typedef IndexSet<3> IndexSetTypeSE;

            // typedef index vector type
            typedef IndexSetTypeSV::ConstIndexVectorReference ConstIndexVectorReferenceSV;
            typedef IndexSetTypeSE::ConstIndexVectorReference ConstIndexVectorReferenceSE;

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
              ConstIndexVectorReferenceSV s_v = index_set_s_v[i];
              ConstIndexVectorReferenceSE s_e = index_set_s_e[i];

              // fetch fine mesh vertices-at-triangle index vectors
              IndexVectorReference s_0 = index_set_out[offset + 4*i + 0];
              IndexVectorReference s_1 = index_set_out[offset + 4*i + 1];
              IndexVectorReference s_2 = index_set_out[offset + 4*i + 2];
              IndexVectorReference s_3 = index_set_out[offset + 4*i + 3];

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
        }; // IndexRefiner<Simplex<2>,2,0>

        /**
         * \brief IndexRefiner implementation for Simplex<2>: Edges-at-Tris
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Simplex<2>, 2, 1>
        {
          enum
          {
            shape_dim = 2,
            cell_dim = 2, // second template parameter
            face_dim = 1  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types;
            typedef IndexSet<2> IndexSetTypeEV; // an edge has 2 vertices
            typedef IndexSet<3> IndexSetTypeSV; // a simplex has 3 vertices
            typedef IndexSet<3> IndexSetTypeSE; // a simplex has 3 edges

            // typedef index vector type
            typedef IndexSetTypeSV::ConstIndexVectorReference ConstIndexVectorReferenceSV;
            typedef IndexSetTypeSE::ConstIndexVectorReference ConstIndexVectorReferenceSE;

            // fetch the index sets
            const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();
            const IndexSetTypeSV& index_set_s_v = index_set_holder_in.get_index_set<2,0>();
            const IndexSetTypeSE& index_set_s_e = index_set_holder_in.get_index_set<2,1>();

            // fetch number of coarse mesh simplices
            const Index num_simps = index_set_s_v.get_num_entities();

            // typedef the sub-index mapping
            typedef Congruency::SubIndexMapping<ShapeType, face_dim, 0> SubIndexMappingType;

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
              ConstIndexVectorReferenceSV s_v = index_set_s_v[i];
              ConstIndexVectorReferenceSE s_e = index_set_s_e[i];

              // create a sub-index mapping object
              SubIndexMappingType sim(s_v, s_e, index_set_e_v);

              // fetch fine mesh edges-at-simplex index vectors
              IndexVectorReference s_0 = index_set_out[offset + 4*i + 0];
              IndexVectorReference s_1 = index_set_out[offset + 4*i + 1];
              IndexVectorReference s_2 = index_set_out[offset + 4*i + 2];
              IndexVectorReference s_3 = index_set_out[offset + 4*i + 3];

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
        }; // IndexRefiner<Simplex<2>,2,1>

        /* ************************************************************************************* */
        /* ************************************************************************************* */
        /* ************************************************************************************* */

        /**
         * \brief IndexRefiner implementation for Simplex<3>: Vertices-at-Edges
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Simplex<3>, 1, 0>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 1, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set type (edges at simplices)
            typedef IndexSet<6> IndexSetTypeSE;

            // typedef index vector type
            typedef IndexSetTypeSE::ConstIndexVectorReference ConstIndexVectorReferenceSE;

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
              ConstIndexVectorReferenceSE s_e = index_set_s_e[i];

              // fetch fine mesh vertices-at-edge index vectors
              IndexVectorReference e_0 = index_set_out[offset + 6*i + 0];
              IndexVectorReference e_1 = index_set_out[offset + 6*i + 1];
              IndexVectorReference e_2 = index_set_out[offset + 6*i + 2];
              IndexVectorReference e_3 = index_set_out[offset + 6*i + 3];
              IndexVectorReference e_4 = index_set_out[offset + 6*i + 4];
              IndexVectorReference e_5 = index_set_out[offset + 6*i + 5];

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
        }; // IndexRefiner<Simplex<3>,1,0>

        /**
         * \brief IndexRefiner implementation for Simplex<3>: Vertices-at-Tris
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Simplex<3>, 2, 0>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 2, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types
            typedef IndexSet<4> IndexSetTypeST;
            typedef IndexSet<6> IndexSetTypeSE;

            // typedef index vector type
            typedef IndexSetTypeST::ConstIndexVectorReference ConstIndexVectorReferenceST;
            typedef IndexSetTypeSE::ConstIndexVectorReference ConstIndexVectorReferenceSE;

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
              ConstIndexVectorReferenceSE s_e = index_set_s_e[i];

              // fetch fine mesh vertices-at-triangle index vectors
              IndexVectorReference t_0 = index_set_out[offset + 16*i + 0];
              IndexVectorReference t_1 = index_set_out[offset + 16*i + 1];
              IndexVectorReference t_2 = index_set_out[offset + 16*i + 2];
              IndexVectorReference t_3 = index_set_out[offset + 16*i + 3];
              IndexVectorReference t_4 = index_set_out[offset + 16*i + 4];
              IndexVectorReference t_5 = index_set_out[offset + 16*i + 5];
              IndexVectorReference t_6 = index_set_out[offset + 16*i + 6];
              IndexVectorReference t_7 = index_set_out[offset + 16*i + 7];
              IndexVectorReference t_8 = index_set_out[offset + 16*i + 8];
              IndexVectorReference t_9 = index_set_out[offset + 16*i + 9];
              IndexVectorReference t_10 = index_set_out[offset + 16*i + 10];
              IndexVectorReference t_11 = index_set_out[offset + 16*i + 11];
              IndexVectorReference t_12 = index_set_out[offset + 16*i + 12];
              IndexVectorReference t_13 = index_set_out[offset + 16*i + 13];
              IndexVectorReference t_14 = index_set_out[offset + 16*i + 14];
              IndexVectorReference t_15 = index_set_out[offset + 16*i + 15];

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
        }; // IndexRefiner<Simplex<3>,2,0>

        /**
         * \brief IndexRefiner implementation for Simplex<3>: Vertices-at-Tetras
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Simplex<3>, 3, 0>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 3, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types
            typedef IndexSet<4> IndexSetTypeST;
            typedef IndexSet<6> IndexSetTypeSE;
            typedef IndexSet<4> IndexSetTypeSV;

            // typedef index vector type
            typedef IndexSetTypeST::ConstIndexVectorReference ConstIndexVectorReferenceST;
            typedef IndexSetTypeSE::ConstIndexVectorReference ConstIndexVectorReferenceSE;
            typedef IndexSetTypeSV::ConstIndexVectorReference ConstIndexVectorReferenceSV;

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
              ConstIndexVectorReferenceSE s_e = index_set_s_e[i];
              ConstIndexVectorReferenceSV s_v = index_set_s_v[i];

              // fetch fine mesh vertices-at-simplex index vectors
              IndexVectorReference s_0 = index_set_out[offset + 12*i + 0];
              IndexVectorReference s_1 = index_set_out[offset + 12*i + 1];
              IndexVectorReference s_2 = index_set_out[offset + 12*i + 2];
              IndexVectorReference s_3 = index_set_out[offset + 12*i + 3];
              IndexVectorReference s_4 = index_set_out[offset + 12*i + 4];
              IndexVectorReference s_5 = index_set_out[offset + 12*i + 5];
              IndexVectorReference s_6 = index_set_out[offset + 12*i + 6];
              IndexVectorReference s_7 = index_set_out[offset + 12*i + 7];
              IndexVectorReference s_8 = index_set_out[offset + 12*i + 8];
              IndexVectorReference s_9 = index_set_out[offset + 12*i + 9];
              IndexVectorReference s_10 = index_set_out[offset + 12*i + 10];
              IndexVectorReference s_11 = index_set_out[offset + 12*i + 11];

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
        }; // IndexRefiner<Simplex<3>,3,0>

        /**
         * \brief IndexRefiner implementation for Simplex<3>: Edges-at-Tris
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Simplex<3>, 2, 1>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 2, // second template parameter
            face_dim = 1  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types;
            typedef IndexSet<3> IndexSetTypeTV;
            typedef IndexSet<4> IndexSetTypeSV;
            typedef IndexSet<4> IndexSetTypeST;

            // typedef index vector type
            typedef IndexSetTypeTV::ConstIndexVectorReference ConstIndexVectorReferenceTV;
            typedef IndexSetTypeSV::ConstIndexVectorReference ConstIndexVectorReferenceSV;
            typedef IndexSetTypeST::ConstIndexVectorReference ConstIndexVectorReferenceST;

            // fetch the index sets
            const IndexSetTypeSV& index_set_s_v = index_set_holder_in.get_index_set<3,0>();
            const IndexSetTypeST& index_set_s_t = index_set_holder_in.get_index_set<3,2>();
            const IndexSetTypeTV& index_set_t_v = index_set_holder_in.get_index_set<2,0>();

            // fetch number of coarse mesh simplices
            const Index num_simps = index_set_s_t.get_num_entities();

            // typedef the sub-index mapping
            typedef Congruency::SubIndexMapping<ShapeType, cell_dim, face_dim> SubIndexMappingType;

            // loop over all coarse mesh simplices
            for(Index i(0); i < num_simps; ++i)
            {
              // fetch coarse mesh vertices-at-simplex and edges-at-simplex index vectors
              ConstIndexVectorReferenceSV s_v = index_set_s_v[i];
              ConstIndexVectorReferenceST s_t = index_set_s_t[i];

              // create a sub-index mapping object
              SubIndexMappingType sim(s_v, s_t, index_set_t_v);

              // fetch fine mesh edges-at-triangle index vectors
              IndexVectorReference t_0 = index_set_out[offset + 16*i + 0];
              IndexVectorReference t_1 = index_set_out[offset + 16*i + 1];
              IndexVectorReference t_2 = index_set_out[offset + 16*i + 2];
              IndexVectorReference t_3 = index_set_out[offset + 16*i + 3];
              IndexVectorReference t_4 = index_set_out[offset + 16*i + 4];
              IndexVectorReference t_5 = index_set_out[offset + 16*i + 5];
              IndexVectorReference t_6 = index_set_out[offset + 16*i + 6];
              IndexVectorReference t_7 = index_set_out[offset + 16*i + 7];
              IndexVectorReference t_8 = index_set_out[offset + 16*i + 8];
              IndexVectorReference t_9 = index_set_out[offset + 16*i + 9];
              IndexVectorReference t_10 = index_set_out[offset + 16*i + 10];
              IndexVectorReference t_11 = index_set_out[offset + 16*i + 11];
              IndexVectorReference t_12 = index_set_out[offset + 16*i + 12];
              IndexVectorReference t_13 = index_set_out[offset + 16*i + 13];
              IndexVectorReference t_14 = index_set_out[offset + 16*i + 14];
              IndexVectorReference t_15 = index_set_out[offset + 16*i + 15];

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
        }; // IndexRefiner<Simplex<3>,2,1>

        /**
         * \brief IndexRefiner implementation for Simplex<3>: Edges-at-Tetras
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Simplex<3>, 3, 1>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 3, // second template parameter
            face_dim = 1  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types
            typedef IndexSet<2> IndexSetTypeEV;
            typedef IndexSet<3> IndexSetTypeTV;
            typedef IndexSet<4> IndexSetTypeSV;
            typedef IndexSet<6> IndexSetTypeSE;
            typedef IndexSet<4> IndexSetTypeST;

            // typedef index vector type
            typedef IndexSetTypeEV::ConstIndexVectorReference ConstIndexVectorReferenceEV;
            typedef IndexSetTypeTV::ConstIndexVectorReference ConstIndexVectorReferenceTV;
            typedef IndexSetTypeSV::ConstIndexVectorReference ConstIndexVectorReferenceSV;
            typedef IndexSetTypeSE::ConstIndexVectorReference ConstIndexVectorReferenceSE;
            typedef IndexSetTypeST::ConstIndexVectorReference ConstIndexVectorReferenceST;

            // fetch the index sets
            const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();
            const IndexSetTypeTV& index_set_t_v = index_set_holder_in.get_index_set<2,0>();
            const IndexSetTypeSV& index_set_s_v = index_set_holder_in.get_index_set<3,0>();
            const IndexSetTypeSE& index_set_s_e = index_set_holder_in.get_index_set<3,1>();
            const IndexSetTypeST& index_set_s_t = index_set_holder_in.get_index_set<3,2>();

            // fetch number of coarse mesh simplices
            const Index num_simps = index_set_s_t.get_num_entities();

            // typedef the sub-index mapping
            typedef Congruency::SubIndexMapping<ShapeType, 2, 1> SubIndexMappingType;
            typedef Congruency::SubIndexMapping<ShapeType, face_dim, 0> EdgeIndexMappingType;

            // loop over all coarse mesh simplices
            for(Index i(0); i < num_simps; ++i)
            {
              // fetch coarse mesh index vectors
              ConstIndexVectorReferenceSV s_v = index_set_s_v[i];
              ConstIndexVectorReferenceSE s_e = index_set_s_e[i];
              ConstIndexVectorReferenceST s_t = index_set_s_t[i];

              // create sub-index mapping objects
              SubIndexMappingType sim(s_v, s_t, index_set_t_v);
              EdgeIndexMappingType edgesim(s_v, s_e, index_set_e_v);

              // fetch fine mesh edges-at-simplex index vectors
              IndexVectorReference s_0 = index_set_out[offset + 12*i + 0];
              IndexVectorReference s_1 = index_set_out[offset + 12*i + 1];
              IndexVectorReference s_2 = index_set_out[offset + 12*i + 2];
              IndexVectorReference s_3 = index_set_out[offset + 12*i + 3];
              IndexVectorReference s_4 = index_set_out[offset + 12*i + 4];
              IndexVectorReference s_5 = index_set_out[offset + 12*i + 5];
              IndexVectorReference s_6 = index_set_out[offset + 12*i + 6];
              IndexVectorReference s_7 = index_set_out[offset + 12*i + 7];
              IndexVectorReference s_8 = index_set_out[offset + 12*i + 8];
              IndexVectorReference s_9 = index_set_out[offset + 12*i + 9];
              IndexVectorReference s_10 = index_set_out[offset + 12*i + 10];
              IndexVectorReference s_11 = index_set_out[offset + 12*i + 11];

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
        }; // IndexRefiner<Simplex<3>,3,1>

        /**
         * \brief IndexRefiner implementation for Simplex<3>: Tris-at-Tetras
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Simplex<3>, 3, 2>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 3, // second template parameter
            face_dim = 2  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types
            typedef IndexSet<3> IndexSetTypeTV;
            typedef IndexSet<4> IndexSetTypeSV;
            typedef IndexSet<4> IndexSetTypeST;

            // typedef index vector type
            typedef IndexSetTypeTV::ConstIndexVectorReference ConstIndexVectorReferenceTV;
            typedef IndexSetTypeSV::ConstIndexVectorReference ConstIndexVectorReferenceSV;
            typedef IndexSetTypeST::ConstIndexVectorReference ConstIndexVectorReferenceST;

            // fetch the index sets
            const IndexSetTypeTV& index_set_t_v = index_set_holder_in.get_index_set<2,0>();
            const IndexSetTypeSV& index_set_s_v = index_set_holder_in.get_index_set<3,0>();
            const IndexSetTypeST& index_set_s_t = index_set_holder_in.get_index_set<3,2>();

            // fetch number of coarse mesh simplices
            const Index num_simps = index_set_s_t.get_num_entities();

            // typedef the sub-index mapping
            typedef Congruency::SubIndexMapping<ShapeType, face_dim, 0> SubIndexMappingType;

            // loop over all coarse mesh simplices
            for(Index i(0); i < num_simps; ++i)
            {
              // fetch coarse mesh vertices-at-simplex and triangle-at-simplex index vectors
              ConstIndexVectorReferenceSV s_v = index_set_s_v[i];
              ConstIndexVectorReferenceST s_t = index_set_s_t[i];

              // create a sub-index mapping object
              SubIndexMappingType sim(s_v, s_t, index_set_t_v);

              // fetch fine mesh triangle-at-simplex index vectors
              IndexVectorReference s_0 = index_set_out[offset + 12*i + 0];
              IndexVectorReference s_1 = index_set_out[offset + 12*i + 1];
              IndexVectorReference s_2 = index_set_out[offset + 12*i + 2];
              IndexVectorReference s_3 = index_set_out[offset + 12*i + 3];
              IndexVectorReference s_4 = index_set_out[offset + 12*i + 4];
              IndexVectorReference s_5 = index_set_out[offset + 12*i + 5];
              IndexVectorReference s_6 = index_set_out[offset + 12*i + 6];
              IndexVectorReference s_7 = index_set_out[offset + 12*i + 7];
              IndexVectorReference s_8 = index_set_out[offset + 12*i + 8];
              IndexVectorReference s_9 = index_set_out[offset + 12*i + 9];
              IndexVectorReference s_10 = index_set_out[offset + 12*i + 10];
              IndexVectorReference s_11 = index_set_out[offset + 12*i + 11];

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
        }; // IndexRefiner<Simplex<3>,3,2>

      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONF_STD_REF_IDX_REF_SIMPLEX_HPP