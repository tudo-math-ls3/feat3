#pragma once
#ifndef KERNEL_GEOMETRY_CONF_STD_REF_IDX_REF_HYPERCUBE_HPP
#define KERNEL_GEOMETRY_CONF_STD_REF_IDX_REF_HYPERCUBE_HPP 1

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
         * \brief IndexRefiner implementation for Hypercube<1>: Vertices-at-Edges
         *
         * \author Peter Zajac
         */
        template<>
        struct IndexRefiner<Shape::Hypercube<1>, 1, 0>
        {
          enum
          {
            shape_dim = 1,
            cell_dim = 1, // second template parameter
            face_dim = 0  // third template parameter
          };

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
        }; // IndexRefiner<Hypercube<1>,1,0>

        /* ************************************************************************************* */
        /* ************************************************************************************* */
        /* ************************************************************************************* */

        /**
         * \brief IndexRefiner implementation for Hypercube<2>: Vertices-at-Edges
         *
         * \author Peter Zajac
         */
        template<>
        struct IndexRefiner<Shape::Hypercube<2>, 1, 0>
        {
          enum
          {
            shape_dim = 2,
            cell_dim = 1, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set type; a quad has 4 edges, so we need an IndexSet<4>
            typedef IndexSet<4> IndexSetTypeQE;

            // typedef index vector type
            typedef IndexSetTypeQE::ConstIndexVectorReference ConstIndexVectorReferenceQE;

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
              ConstIndexVectorReferenceQE q_e = index_set_q_e[i];

              // fetch fine mesh vertices-at-edge index vectors
              IndexVectorReference e_0 = index_set_out[offset + 4*i + 0];
              IndexVectorReference e_1 = index_set_out[offset + 4*i + 1];
              IndexVectorReference e_2 = index_set_out[offset + 4*i + 2];
              IndexVectorReference e_3 = index_set_out[offset + 4*i + 3];

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
        }; // IndexRefiner<Hypercube<2>,1,0>

        /**
         * \brief IndexRefiner implementation for Hypercube<2>: Vertices-at-Quads
         *
         * \author Peter Zajac
         */
        template<>
        struct IndexRefiner<Shape::Hypercube<2>, 2, 0>
        {
          enum
          {
            shape_dim = 2,
            cell_dim = 2, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types; a quad has 4 vertices and 4 edges, so both need an IndexSet<4>
            typedef IndexSet<4> IndexSetTypeQV;
            typedef IndexSet<4> IndexSetTypeQE;

            // typedef index vector type
            typedef IndexSetTypeQV::ConstIndexVectorReference ConstIndexVectorReferenceQV;
            typedef IndexSetTypeQE::ConstIndexVectorReference ConstIndexVectorReferenceQE;

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
              ConstIndexVectorReferenceQV q_v = index_set_q_v[i];
              ConstIndexVectorReferenceQE q_e = index_set_q_e[i];

              // fetch fine mesh vertices-at-quad index vectors
              IndexVectorReference q_0 = index_set_out[offset + 4*i + 0];
              IndexVectorReference q_1 = index_set_out[offset + 4*i + 1];
              IndexVectorReference q_2 = index_set_out[offset + 4*i + 2];
              IndexVectorReference q_3 = index_set_out[offset + 4*i + 3];

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
        }; // IndexRefiner<Hypercube<2>,2,0>

        /**
         * \brief IndexRefiner implementation for Hypercube<2>: Edges-at-Quads
         *
         * \author Peter Zajac
         */
        template<>
        struct IndexRefiner<Shape::Hypercube<2>, 2, 1>
        {
          enum
          {
            shape_dim = 2,
            cell_dim = 2, // second template parameter
            face_dim = 1  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types;
            typedef IndexSet<2> IndexSetTypeEV; // an edge has 2 vertices
            typedef IndexSet<4> IndexSetTypeQV; // a quad has 4 vertices
            typedef IndexSet<4> IndexSetTypeQE; // a quad has 4 edges

            // typedef index vector type
            typedef IndexSetTypeQV::ConstIndexVectorReference ConstIndexVectorReferenceQV;
            typedef IndexSetTypeQE::ConstIndexVectorReference ConstIndexVectorReferenceQE;

            // fetch the index sets
            const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();
            const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
            const IndexSetTypeQE& index_set_q_e = index_set_holder_in.get_index_set<2,1>();

            // fetch number of coarse mesh quads
            const Index num_quads = index_set_q_v.get_num_entities();

            // typedef the sub-index mapping
            typedef Congruency::SubIndexMapping<ShapeType, face_dim, 0> SubIndexMappingType;

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
              ConstIndexVectorReferenceQV q_v = index_set_q_v[i];
              ConstIndexVectorReferenceQE q_e = index_set_q_e[i];

              // create a sub-index mapping object
              SubIndexMappingType sim(q_v, q_e, index_set_e_v);

              // fetch fine mesh edges-at-quad index vectors
              IndexVectorReference q_0 = index_set_out[offset + 4*i + 0];
              IndexVectorReference q_1 = index_set_out[offset + 4*i + 1];
              IndexVectorReference q_2 = index_set_out[offset + 4*i + 2];
              IndexVectorReference q_3 = index_set_out[offset + 4*i + 3];

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
        }; // IndexRefiner<Hypercube<2>,2,1>

        /* ************************************************************************************* */
        /* ************************************************************************************* */
        /* ************************************************************************************* */

        /**
         * \brief IndexRefiner implementation for Hypercube<3>: Vertices-at-Edges
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Hypercube<3>, 1, 0>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 1, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set type
            typedef IndexSet<6> IndexSetTypeCQ;

            // typedef index vector type
            typedef IndexSetTypeCQ::ConstIndexVectorReference ConstIndexVectorReferenceCQ;

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
              ConstIndexVectorReferenceCQ c_q = index_set_c_q[i];

              // fetch fine mesh vertices-at-edge index vectors
              IndexVectorReference e_0 = index_set_out[offset + 6*i + 0];
              IndexVectorReference e_1 = index_set_out[offset + 6*i + 1];
              IndexVectorReference e_2 = index_set_out[offset + 6*i + 2];
              IndexVectorReference e_3 = index_set_out[offset + 6*i + 3];
              IndexVectorReference e_4 = index_set_out[offset + 6*i + 4];
              IndexVectorReference e_5 = index_set_out[offset + 6*i + 5];

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
        }; // IndexRefiner<Hypercube<3>,1,0>

        /**
         * \brief IndexRefiner implementation for Hypercube<3>: Vertices-at-Quads
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Hypercube<3>, 2, 0>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 2, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types
            typedef IndexSet<6> IndexSetTypeCQ;
            typedef IndexSet<12> IndexSetTypeCE;

            // typedef index vector type
            typedef IndexSetTypeCQ::ConstIndexVectorReference ConstIndexVectorReferenceCQ;
            typedef IndexSetTypeCE::ConstIndexVectorReference ConstIndexVectorReferenceCE;

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
              ConstIndexVectorReferenceCQ c_q = index_set_c_q[i];
              ConstIndexVectorReferenceCE c_e = index_set_c_e[i];

              // fetch fine mesh vertices-at-quad index vectors
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
        }; // IndexRefiner<Hypercube<3>,2,0>

        /**
         * \brief IndexRefiner implementation for Hypercube<3>: Vertices-at-Hexas
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Hypercube<3>, 3, 0>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 3, // second template parameter
            face_dim = 0  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types
            typedef IndexSet<6> IndexSetTypeCQ;
            typedef IndexSet<12> IndexSetTypeCE;
            typedef IndexSet<8> IndexSetTypeCV;

            // typedef index vector type
            typedef IndexSetTypeCQ::ConstIndexVectorReference ConstIndexVectorReferenceCQ;
            typedef IndexSetTypeCE::ConstIndexVectorReference ConstIndexVectorReferenceCE;
            typedef IndexSetTypeCV::ConstIndexVectorReference ConstIndexVectorReferenceCV;

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
              ConstIndexVectorReferenceCQ c_q = index_set_c_q[i];
              ConstIndexVectorReferenceCE c_e = index_set_c_e[i];
              ConstIndexVectorReferenceCV c_v = index_set_c_v[i];

              // fetch fine mesh vertices-at-cube index vectors
              IndexVectorReference c_0 = index_set_out[offset + 8*i + 0];
              IndexVectorReference c_1 = index_set_out[offset + 8*i + 1];
              IndexVectorReference c_2 = index_set_out[offset + 8*i + 2];
              IndexVectorReference c_3 = index_set_out[offset + 8*i + 3];
              IndexVectorReference c_4 = index_set_out[offset + 8*i + 4];
              IndexVectorReference c_5 = index_set_out[offset + 8*i + 5];
              IndexVectorReference c_6 = index_set_out[offset + 8*i + 6];
              IndexVectorReference c_7 = index_set_out[offset + 8*i + 7];

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
        }; // IndexRefiner<Hypercube<3>,3,0>

        /**
         * \brief IndexRefiner implementation for Hypercube<3>: Edges-at-Quads
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Hypercube<3>, 2, 1>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 2, // second template parameter
            face_dim = 1  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types;
            typedef IndexSet<4> IndexSetTypeQV;
            typedef IndexSet<8> IndexSetTypeCV;
            typedef IndexSet<6> IndexSetTypeCQ;

            // typedef index vector type
            typedef IndexSetTypeQV::ConstIndexVectorReference ConstIndexVectorReferenceQV;
            typedef IndexSetTypeCV::ConstIndexVectorReference ConstIndexVectorReferenceCV;
            typedef IndexSetTypeCQ::ConstIndexVectorReference ConstIndexVectorReferenceCQ;

            // fetch the index sets
            const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
            const IndexSetTypeCV& index_set_c_v = index_set_holder_in.get_index_set<3,0>();
            const IndexSetTypeCQ& index_set_c_q = index_set_holder_in.get_index_set<3,2>();

            // fetch number of coarse mesh cubes
            const Index num_cubes = index_set_c_q.get_num_entities();

            // typedef the sub-index mapping
            typedef Congruency::SubIndexMapping<ShapeType, cell_dim, face_dim> SubIndexMappingType;

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
              ConstIndexVectorReferenceCV c_v = index_set_c_v[i];
              ConstIndexVectorReferenceCQ c_q = index_set_c_q[i];

              // create a sub-index mapping object
              SubIndexMappingType sim(c_v, c_q, index_set_q_v);

              // fetch fine mesh edges-at-quad index vectors
              IndexVectorReference q_0 = index_set_out[offset + 12*i + 0];
              IndexVectorReference q_1 = index_set_out[offset + 12*i + 1];
              IndexVectorReference q_2 = index_set_out[offset + 12*i + 2];
              IndexVectorReference q_3 = index_set_out[offset + 12*i + 3];
              IndexVectorReference q_4 = index_set_out[offset + 12*i + 4];
              IndexVectorReference q_5 = index_set_out[offset + 12*i + 5];
              IndexVectorReference q_6 = index_set_out[offset + 12*i + 6];
              IndexVectorReference q_7 = index_set_out[offset + 12*i + 7];
              IndexVectorReference q_8 = index_set_out[offset + 12*i + 8];
              IndexVectorReference q_9 = index_set_out[offset + 12*i + 9];
              IndexVectorReference q_10 = index_set_out[offset + 12*i + 10];
              IndexVectorReference q_11 = index_set_out[offset + 12*i + 11];

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
        }; // IndexRefiner<Hypercube<3>,2,1>

        /**
         * \brief IndexRefiner implementation for Hypercube<3>: Edges-at-Hexas
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Hypercube<3>, 3, 1>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 3, // second template parameter
            face_dim = 1  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types
            typedef IndexSet<2> IndexSetTypeEV;
            typedef IndexSet<4> IndexSetTypeQV;
            typedef IndexSet<8> IndexSetTypeCV;
            typedef IndexSet<12> IndexSetTypeCE;
            typedef IndexSet<6> IndexSetTypeCQ;

            // typedef index vector type
            typedef IndexSetTypeEV::ConstIndexVectorReference ConstIndexVectorReferenceEV;
            typedef IndexSetTypeQV::ConstIndexVectorReference ConstIndexVectorReferenceQV;
            typedef IndexSetTypeCV::ConstIndexVectorReference ConstIndexVectorReferenceCV;
            typedef IndexSetTypeCE::ConstIndexVectorReference ConstIndexVectorReferenceCE;
            typedef IndexSetTypeCQ::ConstIndexVectorReference ConstIndexVectorReferenceCQ;

            // fetch the index sets
            const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();
            const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
            const IndexSetTypeCV& index_set_c_v = index_set_holder_in.get_index_set<3,0>();
            const IndexSetTypeCE& index_set_c_e = index_set_holder_in.get_index_set<3,1>();
            const IndexSetTypeCQ& index_set_c_q = index_set_holder_in.get_index_set<3,2>();

            // fetch number of coarse mesh cubes
            const Index num_cubes = index_set_c_q.get_num_entities();

            // typedef the sub-index mapping
            typedef Congruency::SubIndexMapping<ShapeType, 2, 1> SubIndexMappingType;
            typedef Congruency::SubIndexMapping<ShapeType, face_dim, 0> EdgeIndexMappingType;

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
              ConstIndexVectorReferenceCV c_v = index_set_c_v[i];
              ConstIndexVectorReferenceCE c_e = index_set_c_e[i];
              ConstIndexVectorReferenceCQ c_q = index_set_c_q[i];

              // create sub-index mapping objects
              SubIndexMappingType sim(c_v, c_q, index_set_q_v);
              EdgeIndexMappingType edgesim(c_v, c_e, index_set_e_v);

              // fetch fine mesh edges-at-cube index vectors
              IndexVectorReference c_0 = index_set_out[offset + 8*i + 0];
              IndexVectorReference c_1 = index_set_out[offset + 8*i + 1];
              IndexVectorReference c_2 = index_set_out[offset + 8*i + 2];
              IndexVectorReference c_3 = index_set_out[offset + 8*i + 3];
              IndexVectorReference c_4 = index_set_out[offset + 8*i + 4];
              IndexVectorReference c_5 = index_set_out[offset + 8*i + 5];
              IndexVectorReference c_6 = index_set_out[offset + 8*i + 6];
              IndexVectorReference c_7 = index_set_out[offset + 8*i + 7];

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
        }; // IndexRefiner<Hypercube<3>,3,1>

        /**
         * \brief IndexRefiner implementation for Hypercube<3>: Quads-at-Hexas
         *
         * \author Constantin Christof
         */
        template<>
        struct IndexRefiner<Shape::Hypercube<3>, 3, 2>
        {
          enum
          {
            shape_dim = 3,
            cell_dim = 3, // second template parameter
            face_dim = 2  // third template parameter
          };

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
            typedef IndexSetType::IndexVectorReference IndexVectorReference;

            // typedef index set types
            typedef IndexSet<4> IndexSetTypeQV;
            typedef IndexSet<8> IndexSetTypeCV;
            typedef IndexSet<6> IndexSetTypeCQ;

            // typedef index vector type
            typedef IndexSetTypeQV::ConstIndexVectorReference ConstIndexVectorReferenceQV;
            typedef IndexSetTypeCV::ConstIndexVectorReference ConstIndexVectorReferenceCV;
            typedef IndexSetTypeCQ::ConstIndexVectorReference ConstIndexVectorReferenceCQ;

            // fetch the index sets
            const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
            const IndexSetTypeCV& index_set_c_v = index_set_holder_in.get_index_set<3,0>();
            const IndexSetTypeCQ& index_set_c_q = index_set_holder_in.get_index_set<3,2>();

            // fetch number of coarse mesh cubes
            const Index num_cubes = index_set_c_q.get_num_entities();

            // typedef the sub-index mapping
            typedef Congruency::SubIndexMapping<ShapeType, face_dim, 0> SubIndexMappingType;

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
              ConstIndexVectorReferenceCV c_v = index_set_c_v[i];
              ConstIndexVectorReferenceCQ c_q = index_set_c_q[i];

              // create a sub-index mapping object
              SubIndexMappingType sim(c_v, c_q, index_set_q_v);

              // fetch fine mesh quads-at-cube index vectors
              IndexVectorReference c_0 = index_set_out[offset + 8*i + 0];
              IndexVectorReference c_1 = index_set_out[offset + 8*i + 1];
              IndexVectorReference c_2 = index_set_out[offset + 8*i + 2];
              IndexVectorReference c_3 = index_set_out[offset + 8*i + 3];
              IndexVectorReference c_4 = index_set_out[offset + 8*i + 4];
              IndexVectorReference c_5 = index_set_out[offset + 8*i + 5];
              IndexVectorReference c_6 = index_set_out[offset + 8*i + 6];
              IndexVectorReference c_7 = index_set_out[offset + 8*i + 7];

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
        }; // IndexRefiner<Hypercube<3>,3,2>

      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONF_STD_REF_IDX_REF_HYPERCUBE_HPP
