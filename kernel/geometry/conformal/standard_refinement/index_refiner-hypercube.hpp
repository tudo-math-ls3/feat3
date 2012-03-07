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
          typedef IndexSetType::IndexVectorType IdxVecType;
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

            // fetch the edge-vertex index set
            const IndexSetType& index_set_e_v = index_set_holder_in.get_index_set<1,0>();

            // fetch number of coarse mesh edges
            const Index num_edges = index_set_e_v.get_num_entities();

            // loop over all coarse mesh edges
            for(Index i(0); i < num_edges; ++i)
            {
              // fetch coarse mesh edge-vertex index vector
              const IdxVecType& e_v = index_set_e_v[i];

              // fetch fine mesh edge-vertex index vectors
              IdxVecType& idx0 = index_set_out[offset + 2*i + 0];
              IdxVecType& idx1 = index_set_out[offset + 2*i + 1];

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
          typedef IndexSetType::IndexVectorType IdxVecType;
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

            // typedef index set types
            typedef IndexSet<4> IndexSetTypeQE;

            // typedef index vector type
            typedef IndexSetTypeQE::IndexVectorType IndexVectorTypeQE;

            // fetch the quad-vertex index set
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
              // fetch coarse mesh quad-edge index vector
              const IndexVectorTypeQE& q_e = index_set_q_e[i];

              // fetch fine mesh quad-edge index vectors
              IdxVecType& e_0 = index_set_out[offset + 4*i + 0];
              IdxVecType& e_1 = index_set_out[offset + 4*i + 1];
              IdxVecType& e_2 = index_set_out[offset + 4*i + 2];
              IdxVecType& e_3 = index_set_out[offset + 4*i + 3];

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
          typedef IndexSetType::IndexVectorType IdxVecType;
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

            // typedef index set types
            typedef IndexSet<4> IndexSetTypeQV;
            typedef IndexSet<4> IndexSetTypeQE;

            // typedef index vector type
            typedef IndexSetTypeQV::IndexVectorType IndexVectorTypeQV;
            typedef IndexSetTypeQE::IndexVectorType IndexVectorTypeQE;

            // fetch the quad-vertex index set
            const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
            const IndexSetTypeQE& index_set_q_e = index_set_holder_in.get_index_set<2,1>();

            // fetch number of coarse mesh quads
            const Index num_quads = index_set_q_v.get_num_entities();

            // Each coarse mesh quad generates four fine mesh quads (q_i) upon refinement, with
            // each one connecting one of the coarse mesh vertices (v_i), its two neighbour fine
            // mesh vertices generated from edge refinement (e_i) and the vertex generated from
            // quad refinement (q):
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

            // loop over all coarse mesh edges
            for(Index i(0); i < num_quads; ++i)
            {
              // fetch coarse mesh quad-vertex index vector
              const IndexVectorTypeQV& q_v = index_set_q_v[i];
              const IndexVectorTypeQE& q_e = index_set_q_e[i];

              // fetch fine mesh quad-vertex index vectors
              IdxVecType& q_0 = index_set_out[offset + 4*i + 0];
              IdxVecType& q_1 = index_set_out[offset + 4*i + 1];
              IdxVecType& q_2 = index_set_out[offset + 4*i + 2];
              IdxVecType& q_3 = index_set_out[offset + 4*i + 3];

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
          typedef IndexSetType::IndexVectorType IdxVecType;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static Index refine(
            IndexSetType& index_set_out,
            const Index offset,
            const Index index_offsets[],
            const IndexSetHolderType& index_set_holder_in)
          {
            typedef Congruency::SubIndexMapping<ShapeType, face_dim, 0> SubIndexMappingType;

            // fetch vertex index offsets
            const Index ioe = index_offsets[1];
            const Index ioq = index_offsets[2];

            // typedef index set types
            typedef IndexSet<2> IndexSetTypeEV;
            typedef IndexSet<4> IndexSetTypeQV;
            typedef IndexSet<4> IndexSetTypeQE;

            // typedef index vector type
            typedef IndexSetTypeEV::IndexVectorType IndexVectorTypeEV;
            typedef IndexSetTypeQV::IndexVectorType IndexVectorTypeQV;
            typedef IndexSetTypeQE::IndexVectorType IndexVectorTypeQE;

            // fetch the quad-vertex index set
            const IndexSetTypeEV& index_set_e_v = index_set_holder_in.get_index_set<1,0>();
            const IndexSetTypeQV& index_set_q_v = index_set_holder_in.get_index_set<2,0>();
            const IndexSetTypeQE& index_set_q_e = index_set_holder_in.get_index_set<2,1>();

            // fetch number of coarse mesh quads
            const Index num_quads = index_set_q_v.get_num_entities();

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

            // loop over all coarse mesh edges
            for(Index i(0); i < num_quads; ++i)
            {
              // fetch coarse mesh quad-vertex index vector
              const IndexVectorTypeQV& q_v = index_set_q_v[i];
              const IndexVectorTypeQE& q_e = index_set_q_e[i];

              // create a sub-index mapping object
              SubIndexMappingType sim(q_v, q_e, index_set_e_v);

              // fetch fine mesh quad-vertex index vectors
              IdxVecType& q_0 = index_set_out[offset + 4*i + 0];
              IdxVecType& q_1 = index_set_out[offset + 4*i + 1];
              IdxVecType& q_2 = index_set_out[offset + 4*i + 2];
              IdxVecType& q_3 = index_set_out[offset + 4*i + 3];

              // calculate fine quad-vertex indices
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
      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST


#endif // KERNEL_GEOMETRY_CONF_STD_REF_IDX_REF_HYPERCUBE_HPP
