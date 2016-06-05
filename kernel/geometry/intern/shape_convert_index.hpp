#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_SHAPE_CONVERT_INDEX_HPP
#define KERNEL_GEOMETRY_INTERN_SHAPE_CONVERT_INDEX_HPP 1

// includes, FEAT
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/intern/shape_convert_traits.hpp>
#include <kernel/geometry/intern/entity_counter.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Shape_,
        int cell_dim_>
      struct ShapeConvertIndex;

      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* ***************************** Hypercube --> Simplex *********************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */

      /**
       * \brief ShapeConvertIndex implementation for Hypercube<1> --> Simplex<1>
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertIndex<Shape::Simplex<1>, 1>
      {
        static constexpr int shape_dim = 1;
        static constexpr int cell_dim = 1;
        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::Hypercube<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<2> IndexSetTypeIn;

          // typedef index vector type
          typedef IndexSetTypeIn::ConstIndexVectorReference ConstIndexVectorReference;

          // fetch the vertices-at-edge index set
          const IndexSetTypeIn& index_set_in = index_set_holder_in.get_index_set<1,0>();

          // fetch number of coarse mesh edges
          const Index num_edges = index_set_in.get_num_entities();

          // loop over all coarse mesh quads
          for(Index i(0); i < num_edges; ++i)
          {
            // fetch coarse mesh vertices-at-quad index vector
            ConstIndexVectorReference v_e = index_set_in[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexVectorReference e_0 = index_set_out[offset + i];

            e_0[0] = iov + v_e[0];
            e_0[1] = iov + v_e[1];
          }

          // return inner edge count
          return num_edges;
        }
      };

      /**
       * \brief ShapeConvertIndex implementation for Hypercube<2> --> Simplex<1>
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertIndex<Shape::Simplex<2>, 1>
      {
        static constexpr int shape_dim = 2;
        static constexpr int cell_dim = 1;
        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::Hypercube<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioq = index_offsets[2];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<4> IndexSetTypeIn;

          // typedef index vector type
          typedef IndexSetTypeIn::ConstIndexVectorReference ConstIndexVectorReference;

          // fetch the vertices-at-quad index set
          const IndexSetTypeIn& index_set_in = index_set_holder_in.get_index_set<2,0>();

          // fetch number of coarse mesh quads
          const Index num_quads = index_set_in.get_num_entities();

          // Each quad generates four inner edges (e_i) upon conversion,
          // with each one connecting one of the quad's corner vertices (v_i)
          // with the vertex generated from the quad midpoint (X).
          //
          // v_2-------------------------v_3
          //  |\__                     __/|
          //  |   \                   /   |
          //  |    e_2__         __e_3    |
          //  |         \__   __/         |
          //  |            \ /            |
          //  |             X             |
          //  |          __/ \__          |
          //  |       __/       \__       |
          //  |    e_0             e_1    |
          //  | __/                   \__ |
          //  |/                         \|
          // v_0-------------------------v_1


          // loop over all coarse mesh quads
          for(Index i(0); i < num_quads; ++i)
          {
            // fetch coarse mesh vertices-at-quad index vector
            ConstIndexVectorReference v_q = index_set_in[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexVectorReference e_0 = index_set_out[offset + 4*i + 0];
            IndexVectorReference e_1 = index_set_out[offset + 4*i + 1];
            IndexVectorReference e_2 = index_set_out[offset + 4*i + 2];
            IndexVectorReference e_3 = index_set_out[offset + 4*i + 3];

            e_0[0] = iov + v_q[0]; // v_0
            e_0[1] = ioq + i; // X

            e_1[0] = iov + v_q[1]; // v_1
            e_1[1] = ioq + i; // X

            e_2[0] = iov + v_q[2]; // v_2
            e_2[1] = ioq + i; // X

            e_3[0] = iov + v_q[3]; // v_3
            e_3[1] = ioq + i; // X
          }

          // return inner edge count
          return 4*num_quads;
        }
      };

      /**
       * \brief ShapeConvertIndex implementation for Hypercube<3> --> Simplex<1>
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertIndex<Shape::Simplex<3>, 1>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 1;
        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::Hypercube<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioq = index_offsets[2];
          const Index ioc = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<8> IndexSetTypeIn_v_c;
          typedef IndexSet<6> IndexSetTypeIn_q_c;

          // typedef index vector type
          typedef IndexSetTypeIn_v_c::ConstIndexVectorReference ConstIndexVectorReference_v_c;
          typedef IndexSetTypeIn_q_c::ConstIndexVectorReference ConstIndexVectorReference_q_c;

          // fetch the vertices-at-hexa and the quad-at-hexa index set
          const IndexSetTypeIn_v_c& index_set_in_v_c = index_set_holder_in.get_index_set<3,0>();
          const IndexSetTypeIn_q_c& index_set_in_q_c = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh hexas
          const Index num_hexas = index_set_in_v_c.get_num_entities();

          // loop over all coarse mesh hexas
          for(Index i(0); i < num_hexas; ++i)
          {
            // fetch coarse mesh vertices-at-hexa index vector
            ConstIndexVectorReference_v_c v_c = index_set_in_v_c[i];
            ConstIndexVectorReference_q_c q_c = index_set_in_q_c[i];

            // fetch new mesh vertices-at-edge index vectors
            IndexVectorReference e_0 = index_set_out[offset + 14*i + 0];
            IndexVectorReference e_1 = index_set_out[offset + 14*i + 1];
            IndexVectorReference e_2 = index_set_out[offset + 14*i + 2];
            IndexVectorReference e_3 = index_set_out[offset + 14*i + 3];
            IndexVectorReference e_4 = index_set_out[offset + 14*i + 4];
            IndexVectorReference e_5 = index_set_out[offset + 14*i + 5];
            IndexVectorReference e_6 = index_set_out[offset + 14*i + 6];
            IndexVectorReference e_7 = index_set_out[offset + 14*i + 7];
            IndexVectorReference e_8 = index_set_out[offset + 14*i + 8];
            IndexVectorReference e_9 = index_set_out[offset + 14*i + 9];
            IndexVectorReference e_10 = index_set_out[offset + 14*i + 10];
            IndexVectorReference e_11 = index_set_out[offset + 14*i + 11];
            IndexVectorReference e_12 = index_set_out[offset + 14*i + 12];
            IndexVectorReference e_13 = index_set_out[offset + 14*i + 13];

            // corner-midpoint edges
            e_0[0] = iov + v_c[0]; // v_0
            e_0[1] = ioc + i; // X

            e_1[0] = iov + v_c[1]; // v_1
            e_1[1] = ioc + i; // X

            e_2[0] = iov + v_c[2]; // v_2
            e_2[1] = ioc + i; // X

            e_3[0] = iov + v_c[3]; // v_3
            e_3[1] = ioc + i; // X

            e_4[0] = iov + v_c[4]; // v_4
            e_4[1] = ioc + i; // X

            e_5[0] = iov + v_c[5]; // v_5
            e_5[1] = ioc + i; // X

            e_6[0] = iov + v_c[6]; // v_6
            e_6[1] = ioc + i; // X

            e_7[0] = iov + v_c[7]; // v_7
            e_7[1] = ioc + i; // X


            // face-midpoint edges
            e_8[0] = ioq + q_c[0]; // q_0
            e_8[1] = ioc + i; // X

            e_9[0] = ioq + q_c[1]; // q_1
            e_9[1] = ioc + i; // X

            e_10[0] = ioq + q_c[2]; // q_2
            e_10[1] = ioc + i; // X

            e_11[0] = ioq + q_c[3]; // q_3
            e_11[1] = ioc + i; // X

            e_12[0] = ioq + q_c[4]; // q_4
            e_12[1] = ioc + i; // X

            e_13[0] = ioq + q_c[5]; // q_5
            e_13[1] = ioc + i; // X

          }

          // return inner edge count
          return 14*num_hexas;
        }
      };

      /**
       * \brief ShapeConvertIndex implementation for Hypercube<2> --> Simplex<2>
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertIndex<Shape::Simplex<2>, 2>
      {
        static constexpr int shape_dim = 2;
        static constexpr int cell_dim = 2;
        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::Hypercube<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioq = index_offsets[2];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<4> IndexSetTypeIn;

          // typedef index vector type
          typedef IndexSetTypeIn::ConstIndexVectorReference ConstIndexVectorReference;

          // fetch the vertices-at-quad index set
          const IndexSetTypeIn& index_set_in = index_set_holder_in.get_index_set<2,0>();

          // fetch number of coarse mesh quads
          const Index num_quads = index_set_in.get_num_entities();

          // Each quad generates four inner triangles (t_i) upon conversion,
          // with each one connecting two of the quad's corner vertices (v_i)
          // with the vertex generated from the quad midpoint (X).
          //
          // v_2-------------------------v_3
          //  |\__                     __/|
          //  |   \__      t_1      __/   |
          //  |      \__         __/      |
          //  |         \__   __/         |
          //  |            \ /            |
          //  |    t_2      X      t_3    |
          //  |          __/ \__          |
          //  |       __/       \__       |
          //  |    __/             \__    |
          //  | __/        t_0        \__ |
          //  |/                         \|
          // v_0-------------------------v_1


          // loop over all coarse mesh quads
          for(Index i(0); i < num_quads; ++i)
          {
            // fetch coarse mesh vertices-at-quad index vector
            ConstIndexVectorReference v_q = index_set_in[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexVectorReference t_0 = index_set_out[offset + 4*i + 0];
            IndexVectorReference t_1 = index_set_out[offset + 4*i + 1];
            IndexVectorReference t_2 = index_set_out[offset + 4*i + 2];
            IndexVectorReference t_3 = index_set_out[offset + 4*i + 3];

            t_0[0] = ioq + i; // X
            t_0[1] = iov + v_q[0]; // v_0
            t_0[2] = iov + v_q[1]; // v_1

            t_1[0] = ioq + i; // X
            t_1[1] = iov + v_q[3]; // v_3
            t_1[2] = iov + v_q[2]; // v_2

            t_2[0] = ioq + i; // X
            t_2[1] = iov + v_q[2]; // v_2
            t_2[2] = iov + v_q[0]; // v_0

            t_3[0] = ioq + i; // X
            t_3[1] = iov + v_q[1]; // v_1
            t_3[2] = iov + v_q[3]; // v_3
          }

          // return triangle count
          return 4*num_quads;
        }
      };


      /**
       * \brief ShapeConvertIndex implementation for Hypercube<3> --> Simplex<2>
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertIndex<Shape::Simplex<3>, 2>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 2;
        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::Hypercube<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioq = index_offsets[2];
          const Index ioc = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<8> IndexSetTypeIn_v_c;
          typedef IndexSet<6> IndexSetTypeIn_q_c;

          // typedef index vector type
          typedef IndexSetTypeIn_v_c::ConstIndexVectorReference ConstIndexVectorReference_v_c;
          typedef IndexSetTypeIn_q_c::ConstIndexVectorReference ConstIndexVectorReference_q_c;

          // fetch the vertices-at-hexa and quad-at-hexa index set
          const IndexSetTypeIn_v_c& index_set_in_v_c = index_set_holder_in.get_index_set<3,0>();
          const IndexSetTypeIn_q_c& index_set_in_q_c = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh hexas
          const Index num_cubes = index_set_in_v_c.get_num_entities();

          // loop over all coarse mesh hexas
          for(Index i(0); i < num_cubes; ++i)
          {
            // fetch coarse mesh vertices-at-hexas and quad-at-hexa index vector
            ConstIndexVectorReference_v_c v_c = index_set_in_v_c[i];
            ConstIndexVectorReference_q_c q_c = index_set_in_q_c[i];

            // fetch fine mesh vertices-at-triangle index vectors
            IndexVectorReference t_0 = index_set_out[offset + 36*i + 0];
            IndexVectorReference t_1 = index_set_out[offset + 36*i + 1];
            IndexVectorReference t_2 = index_set_out[offset + 36*i + 2];
            IndexVectorReference t_3 = index_set_out[offset + 36*i + 3];
            IndexVectorReference t_4 = index_set_out[offset + 36*i + 4];
            IndexVectorReference t_5 = index_set_out[offset + 36*i + 5];
            IndexVectorReference t_6 = index_set_out[offset + 36*i + 6];
            IndexVectorReference t_7 = index_set_out[offset + 36*i + 7];
            IndexVectorReference t_8 = index_set_out[offset + 36*i + 8];
            IndexVectorReference t_9 = index_set_out[offset + 36*i + 9];
            IndexVectorReference t_10 = index_set_out[offset + 36*i + 10];
            IndexVectorReference t_11 = index_set_out[offset + 36*i + 11];
            IndexVectorReference t_12 = index_set_out[offset + 36*i + 12];
            IndexVectorReference t_13 = index_set_out[offset + 36*i + 13];
            IndexVectorReference t_14 = index_set_out[offset + 36*i + 14];
            IndexVectorReference t_15 = index_set_out[offset + 36*i + 15];
            IndexVectorReference t_16 = index_set_out[offset + 36*i + 16];
            IndexVectorReference t_17 = index_set_out[offset + 36*i + 17];
            IndexVectorReference t_18 = index_set_out[offset + 36*i + 18];
            IndexVectorReference t_19 = index_set_out[offset + 36*i + 19];
            IndexVectorReference t_20 = index_set_out[offset + 36*i + 20];
            IndexVectorReference t_21 = index_set_out[offset + 36*i + 21];
            IndexVectorReference t_22 = index_set_out[offset + 36*i + 22];
            IndexVectorReference t_23 = index_set_out[offset + 36*i + 23];
            IndexVectorReference t_24 = index_set_out[offset + 36*i + 24];
            IndexVectorReference t_25 = index_set_out[offset + 36*i + 25];
            IndexVectorReference t_26 = index_set_out[offset + 36*i + 26];
            IndexVectorReference t_27 = index_set_out[offset + 36*i + 27];
            IndexVectorReference t_28 = index_set_out[offset + 36*i + 28];
            IndexVectorReference t_29 = index_set_out[offset + 36*i + 29];
            IndexVectorReference t_30 = index_set_out[offset + 36*i + 30];
            IndexVectorReference t_31 = index_set_out[offset + 36*i + 31];
            IndexVectorReference t_32 = index_set_out[offset + 36*i + 32];
            IndexVectorReference t_33 = index_set_out[offset + 36*i + 33];
            IndexVectorReference t_34 = index_set_out[offset + 36*i + 34];
            IndexVectorReference t_35 = index_set_out[offset + 36*i + 35];

            // coarse edge triangles
            t_0[0] = ioc + i; // X
            t_0[1] = iov + v_c[0]; // v_0
            t_0[2] = iov + v_c[1]; // v_1

            t_1[0] = ioc + i; // X
            t_1[1] = iov + v_c[3]; // v_3
            t_1[2] = iov + v_c[2]; // v_2

            t_2[0] = ioc + i; // X
            t_2[1] = iov + v_c[4]; // v_4
            t_2[2] = iov + v_c[5]; // v_5

            t_3[0] = ioc + i; // X
            t_3[1] = iov + v_c[7]; // v_7
            t_3[2] = iov + v_c[6]; // v_6

            t_4[0] = ioc + i; // X
            t_4[1] = iov + v_c[2]; // v_2
            t_4[2] = iov + v_c[0]; // v_0

            t_5[0] = ioc + i; // X
            t_5[1] = iov + v_c[1]; // v_1
            t_5[2] = iov + v_c[3]; // v_3

            t_6[0] = ioc + i; // X
            t_6[1] = iov + v_c[6]; // v_6
            t_6[2] = iov + v_c[4]; // v_4

            t_7[0] = ioc + i; // X
            t_7[1] = iov + v_c[5]; // v_5
            t_7[2] = iov + v_c[7]; // v_7

            t_8[0] = ioc + i; // X
            t_8[1] = iov + v_c[4]; // v_4
            t_8[2] = iov + v_c[0]; // v_0

            t_9[0] = ioc + i; // X
            t_9[1] = iov + v_c[1]; // v_1
            t_9[2] = iov + v_c[5]; // v_5

            t_10[0] = ioc + i; // X
            t_10[1] = iov + v_c[6]; // v_6
            t_10[2] = iov + v_c[2]; // v_2

            t_11[0] = ioc + i; // X
            t_11[1] = iov + v_c[3]; // v_3
            t_11[2] = iov + v_c[7]; // v_7

            // face triangles
            // face 0
            t_12[0] = ioc + i; // X
            t_12[1] = ioq + q_c[0]; // q_0
            t_12[2] = iov + v_c[0]; // v_0

            t_13[0] = ioc + i; // X
            t_13[1] = ioq + q_c[0]; // q_0
            t_13[2] = iov + v_c[1]; // v_1

            t_14[0] = ioc + i; // X
            t_14[1] = ioq + q_c[0]; // q_0
            t_14[2] = iov + v_c[2]; // v_2

            t_15[0] = ioc + i; // X
            t_15[1] = ioq + q_c[0]; // q_0
            t_15[2] = iov + v_c[3]; // v_3

            // face 1
            t_16[0] = ioc + i; // X
            t_16[1] = ioq + q_c[1]; // q_1
            t_16[2] = iov + v_c[4]; // v_4

            t_17[0] = ioc + i; // X
            t_17[1] = ioq + q_c[1]; // q_1
            t_17[2] = iov + v_c[5]; // v_5

            t_18[0] = ioc + i; // X
            t_18[1] = ioq + q_c[1]; // q_1
            t_18[2] = iov + v_c[6]; // v_6

            t_19[0] = ioc + i; // X
            t_19[1] = ioq + q_c[1]; // q_1
            t_19[2] = iov + v_c[7]; // v_7

            // face 2
            t_20[0] = ioc + i; // X
            t_20[1] = ioq + q_c[2]; // q_2
            t_20[2] = iov + v_c[0]; // v_0

            t_21[0] = ioc + i; // X
            t_21[1] = ioq + q_c[2]; // q_2
            t_21[2] = iov + v_c[1]; // v_1

            t_22[0] = ioc + i; // X
            t_22[1] = ioq + q_c[2]; // q_2
            t_22[2] = iov + v_c[4]; // v_4

            t_23[0] = ioc + i; // X
            t_23[1] = ioq + q_c[2]; // q_2
            t_23[2] = iov + v_c[5]; // v_5

            // face 3
            t_24[0] = ioc + i; // X
            t_24[1] = ioq + q_c[3]; // q_3
            t_24[2] = iov + v_c[2]; // v_2

            t_25[0] = ioc + i; // X
            t_25[1] = ioq + q_c[3]; // q_3
            t_25[2] = iov + v_c[3]; // v_3

            t_26[0] = ioc + i; // X
            t_26[1] = ioq + q_c[3]; // q_3
            t_26[2] = iov + v_c[6]; // v_6

            t_27[0] = ioc + i; // X
            t_27[1] = ioq + q_c[3]; // q_3
            t_27[2] = iov + v_c[7]; // v_7

            // face 4
            t_28[0] = ioc + i; // X
            t_28[1] = ioq + q_c[4]; // q_4
            t_28[2] = iov + v_c[0]; // v_0

            t_29[0] = ioc + i; // X
            t_29[1] = ioq + q_c[4]; // q_4
            t_29[2] = iov + v_c[2]; // v_2

            t_30[0] = ioc + i; // X
            t_30[1] = ioq + q_c[4]; // q_4
            t_30[2] = iov + v_c[4]; // v_4

            t_31[0] = ioc + i; // X
            t_31[1] = ioq + q_c[4]; // q_4
            t_31[2] = iov + v_c[6]; // v_6

            // face 5
            t_32[0] = ioc + i; // X
            t_32[1] = ioq + q_c[5]; // q_5
            t_32[2] = iov + v_c[1]; // v_1

            t_33[0] = ioc + i; // X
            t_33[1] = ioq + q_c[5]; // q_5
            t_33[2] = iov + v_c[3]; // v_3

            t_34[0] = ioc + i; // X
            t_34[1] = ioq + q_c[5]; // q_5
            t_34[2] = iov + v_c[5]; // v_5

            t_35[0] = ioc + i; // X
            t_35[1] = ioq + q_c[5]; // q_5
            t_35[2] = iov + v_c[7]; // v_7

          }

          // return triangle count
          return 36*num_cubes;
        }
      };


      /**
       * \brief ShapeConvertIndex implementation for Hypercube<3> --> Simplex<3>
       *
       * \author Stefan Wahlers, Constantin Christof
       */

      template<>
      struct ShapeConvertIndex<Shape::Simplex<3>, 3>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 3;
        typedef Shape::Simplex<shape_dim> ShapeType;
        typedef Shape::Hypercube<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioq = index_offsets[2];
          const Index ioc = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<8> IndexSetTypeIn_v_c;
          typedef IndexSet<6> IndexSetTypeIn_q_c;

          // typedef index vector type
          typedef IndexSetTypeIn_v_c::ConstIndexVectorReference ConstIndexVectorReference_v_c;
          typedef IndexSetTypeIn_q_c::ConstIndexVectorReference ConstIndexVectorReference_q_c;

          // fetch the vertices-at-hexas and quad-at-hexas index set
          const IndexSetTypeIn_v_c& index_set_in_v_c = index_set_holder_in.get_index_set<3,0>();
          const IndexSetTypeIn_q_c& index_set_in_q_c = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh hexas
          const Index num_cubes = index_set_in_v_c.get_num_entities();

          // loop over all coarse mesh hexas
          for(Index i(0); i < num_cubes; ++i)
          {
            // fetch coarse mesh vertices-at-hexas and quad-at-hexas index vector
            ConstIndexVectorReference_v_c v_c = index_set_in_v_c[i];
            ConstIndexVectorReference_q_c q_c = index_set_in_q_c[i];

            // fetch fine mesh vertices-at-tetra index vectors
            IndexVectorReference t_0 = index_set_out[offset + 24*i + 0];
            IndexVectorReference t_1 = index_set_out[offset + 24*i + 1];
            IndexVectorReference t_2 = index_set_out[offset + 24*i + 2];
            IndexVectorReference t_3 = index_set_out[offset + 24*i + 3];
            IndexVectorReference t_4 = index_set_out[offset + 24*i + 4];
            IndexVectorReference t_5 = index_set_out[offset + 24*i + 5];
            IndexVectorReference t_6 = index_set_out[offset + 24*i + 6];
            IndexVectorReference t_7 = index_set_out[offset + 24*i + 7];
            IndexVectorReference t_8 = index_set_out[offset + 24*i + 8];
            IndexVectorReference t_9 = index_set_out[offset + 24*i + 9];
            IndexVectorReference t_10 = index_set_out[offset + 24*i + 10];
            IndexVectorReference t_11 = index_set_out[offset + 24*i + 11];
            IndexVectorReference t_12 = index_set_out[offset + 24*i + 12];
            IndexVectorReference t_13 = index_set_out[offset + 24*i + 13];
            IndexVectorReference t_14 = index_set_out[offset + 24*i + 14];
            IndexVectorReference t_15 = index_set_out[offset + 24*i + 15];
            IndexVectorReference t_16 = index_set_out[offset + 24*i + 16];
            IndexVectorReference t_17 = index_set_out[offset + 24*i + 17];
            IndexVectorReference t_18 = index_set_out[offset + 24*i + 18];
            IndexVectorReference t_19 = index_set_out[offset + 24*i + 19];
            IndexVectorReference t_20 = index_set_out[offset + 24*i + 20];
            IndexVectorReference t_21 = index_set_out[offset + 24*i + 21];
            IndexVectorReference t_22 = index_set_out[offset + 24*i + 22];
            IndexVectorReference t_23 = index_set_out[offset + 24*i + 23];

            // end up with X and use convention for tetrahedron index numbering (see wiki)
            // face 0
            t_0[0] = iov + v_c[0]; // v_0
            t_0[1] = iov + v_c[1]; // v_1
            t_0[2] = ioq + q_c[0]; // q_0
            t_0[3] = ioc + i; // X

            t_1[0] = iov + v_c[3]; // v_3
            t_1[1] = iov + v_c[2]; // v_2
            t_1[2] = ioq + q_c[0]; // q_0
            t_1[3] = ioc + i; // X

            t_2[0] = iov + v_c[2]; // v_2
            t_2[1] = iov + v_c[0]; // v_0
            t_2[2] = ioq + q_c[0]; // q_0
            t_2[3] = ioc + i; // X

            t_3[0] = iov + v_c[1]; // v_1
            t_3[1] = iov + v_c[3]; // v_3
            t_3[2] = ioq + q_c[0]; // q_0
            t_3[3] = ioc + i; // X

            // face 1
            t_4[0] = iov + v_c[5]; // v_5
            t_4[1] = iov + v_c[4]; // v_4
            t_4[2] = ioq + q_c[1]; // q_1
            t_4[3] = ioc + i; // X

            t_5[0] = iov + v_c[6]; // v_6
            t_5[1] = iov + v_c[7]; // v_7
            t_5[2] = ioq + q_c[1]; // q_1
            t_5[3] = ioc + i; // X

            t_6[0] = iov + v_c[4]; // v_4
            t_6[1] = iov + v_c[6]; // v_6
            t_6[2] = ioq + q_c[1]; // q_1
            t_6[3] = ioc + i; // X

            t_7[0] = iov + v_c[7]; // v_7
            t_7[1] = iov + v_c[5]; // v_5
            t_7[2] = ioq + q_c[1]; // q_1
            t_7[3] = ioc + i; // X

            // face 2
            t_8[0] = iov + v_c[1]; // v_1
            t_8[1] = iov + v_c[0]; // v_0
            t_8[2] = ioq + q_c[2]; // q_2
            t_8[3] = ioc + i; // X

            t_9[0] = iov + v_c[4]; // v_4
            t_9[1] = iov + v_c[5]; // v_5
            t_9[2] = ioq + q_c[2]; // q_2
            t_9[3] = ioc + i; // X

            t_10[0] = iov + v_c[0]; // v_0
            t_10[1] = iov + v_c[4]; // v_4
            t_10[2] = ioq + q_c[2]; // q_2
            t_10[3] = ioc + i; // X

            t_11[0] = iov + v_c[5]; // v_5
            t_11[1] = iov + v_c[1]; // v_1
            t_11[2] = ioq + q_c[2]; // q_2
            t_11[3] = ioc + i; // X

            // face 3
            t_12[0] = iov + v_c[2]; // v_2
            t_12[1] = iov + v_c[3]; // v_3
            t_12[2] = ioq + q_c[3]; // q_3
            t_12[3] = ioc + i; // X

            t_13[0] = iov + v_c[7]; // v_7
            t_13[1] = iov + v_c[6]; // v_6
            t_13[2] = ioq + q_c[3]; // q_3
            t_13[3] = ioc + i; // X

            t_14[0] = iov + v_c[6]; // v_6
            t_14[1] = iov + v_c[2]; // v_2
            t_14[2] = ioq + q_c[3]; // q_3
            t_14[3] = ioc + i; // X

            t_15[0] = iov + v_c[3]; // v_3
            t_15[1] = iov + v_c[7]; // v_7
            t_15[2] = ioq + q_c[3]; // q_3
            t_15[3] = ioc + i; // X

            // face 4
            t_16[0] = iov + v_c[0]; // v_0
            t_16[1] = iov + v_c[2]; // v_2
            t_16[2] = ioq + q_c[4]; // q_4
            t_16[3] = ioc + i; // X

            t_17[0] = iov + v_c[6]; // v_6
            t_17[1] = iov + v_c[4]; // v_4
            t_17[2] = ioq + q_c[4]; // q_4
            t_17[3] = ioc + i; // X

            t_18[0] = iov + v_c[4]; // v_4
            t_18[1] = iov + v_c[0]; // v_0
            t_18[2] = ioq + q_c[4]; // q_4
            t_18[3] = ioc + i; // X

            t_19[0] = iov + v_c[2]; // v_2
            t_19[1] = iov + v_c[6]; // v_6
            t_19[2] = ioq + q_c[4]; // q_4
            t_19[3] = ioc + i; // X

            // face 5
            t_20[0] = iov + v_c[3]; // v_3
            t_20[1] = iov + v_c[1]; // v_1
            t_20[2] = ioq + q_c[5]; // q_5
            t_20[3] = ioc + i; // X

            t_21[0] = iov + v_c[5]; // v_5
            t_21[1] = iov + v_c[7]; // v_7
            t_21[2] = ioq + q_c[5]; // q_5
            t_21[3] = ioc + i; // X

            t_22[0] = iov + v_c[1]; // v_1
            t_22[1] = iov + v_c[5]; // v_5
            t_22[2] = ioq + q_c[5]; // q_5
            t_22[3] = ioc + i; // X

            t_23[0] = iov + v_c[7]; // v_7
            t_23[1] = iov + v_c[3]; // v_3
            t_23[2] = ioq + q_c[5]; // q_5
            t_23[3] = ioc + i; // X
          }

          // return triangle count
          return 24*num_cubes;
        }
      };

      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* ***************************** Simplex --> Hypercube *********************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */

      /**
       * \brief ShapeConvertIndex implementation for Simplex<1> --> Hypercube<1>
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertIndex<Shape::Hypercube<1>, 1>
      {
        static constexpr int shape_dim = 1;
        static constexpr int cell_dim = 1;
        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::Simplex<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

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

          // typedef index set types
          typedef IndexSet<2> IndexSetTypeIn;

          // typedef index vector type
          typedef IndexSetTypeIn::ConstIndexVectorReference ConstIndexVectorReference;

          // fetch the vertices-at-edge index set
          const IndexSetTypeIn& index_set_in = index_set_holder_in.get_index_set<1,0>();

          // fetch number of coarse mesh edges
          const Index num_edges = index_set_in.get_num_entities();

          // loop over all coarse mesh edges
          for(Index i(0); i < num_edges; ++i)
          {
            // fetch coarse mesh vertices-at-edge index vector
            ConstIndexVectorReference v_e = index_set_in[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexVectorReference e_0 = index_set_out[offset + 2*i+0];
            IndexVectorReference e_1 = index_set_out[offset + 2*i+1];

            e_0[0] = iov + v_e[0];
            e_0[1] = ioe + i;

            e_1[0] = ioe + i;
            e_1[1] = iov + v_e[1];
          }

          // return inner edge count
          return 2*num_edges;
        }
      };

      /**
       * \brief ShapeConvertIndex implementation for Simplex<2> --> Hypercube<1>
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertIndex<Shape::Hypercube<2>, 1>
      {
        static constexpr int shape_dim = 2;
        static constexpr int cell_dim = 1;
        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::Simplex<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index ioe = index_offsets[1];
          const Index iot = index_offsets[2];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<3> IndexSetTypeIn;

          // typedef index vector type
          typedef IndexSetTypeIn::ConstIndexVectorReference ConstIndexVectorReference;

          // fetch the edge-at-triangle index set
          const IndexSetTypeIn& index_set_in = index_set_holder_in.get_index_set<2,1>();

          // fetch number of coarse mesh triangle
          const Index num_trias = index_set_in.get_num_entities();

          // loop over all coarse mesh triangle
          for(Index i(0); i < num_trias; ++i)
          {
            // fetch coarse mesh edge-at-triangle index vector
            ConstIndexVectorReference e_t = index_set_in[i];

            // fetch fine mesh vertices-at-edge index vectors
            IndexVectorReference e_0 = index_set_out[offset + 3*i + 0];
            IndexVectorReference e_1 = index_set_out[offset + 3*i + 1];
            IndexVectorReference e_2 = index_set_out[offset + 3*i + 2];

            // end with X
            e_0[0] = ioe + e_t[0]; // v_0
            e_0[1] = iot + i; // X

            e_1[0] = ioe + e_t[1]; // v_1
            e_1[1] = iot + i; // X

            e_2[0] = ioe + e_t[2]; // v_2
            e_2[1] = iot + i; // X
          }

          // return inner edge count
          return 3*num_trias;
        }
      };

      /**
       * \brief ShapeConvertIndex implementation for Simplex<3> --> Hypercube<1>
       *
       * \author: Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertIndex<Shape::Hypercube<3>, 1>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 1;
        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::Simplex<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iot = index_offsets[2];
          const Index iotet = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<4> IndexSetTypeIn_t_tet;

          // typedef index vector type
          typedef IndexSetTypeIn_t_tet::ConstIndexVectorReference ConstIndexVectorReference_t_tet;

          // fetch the triangle-at-tetra
          const IndexSetTypeIn_t_tet& index_set_in_t_tet = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh tetras
          const Index num_tetras = index_set_in_t_tet.get_num_entities();

          // loop over all coarse mesh tetras
          for(Index i(0); i < num_tetras; ++i)
          {
            // fetch coarse mesh triangle-at-tetra index vector
            ConstIndexVectorReference_t_tet t_tet = index_set_in_t_tet[i];

            // fetch new mesh vertices-at-edge index vectors
            IndexVectorReference e_0 = index_set_out[offset + 4*i + 0];
            IndexVectorReference e_1 = index_set_out[offset + 4*i + 1];
            IndexVectorReference e_2 = index_set_out[offset + 4*i + 2];
            IndexVectorReference e_3 = index_set_out[offset + 4*i + 3];


            // end with X
            e_0[0] = iot + t_tet[0]; // v_0
            e_0[1] = iotet + i; // X

            e_1[0] = iot + t_tet[1]; // v_1
            e_1[1] = iotet + i; // X

            e_2[0] = iot + t_tet[2]; // v_2
            e_2[1] = iotet + i; // X

            e_3[0] = iot + t_tet[3]; // v_3
            e_3[1] = iotet + i; // X

          }

          // return inner edge count
          return 4*num_tetras;
        }
      };

      /**
       * \brief ShapeConvertIndex implementation for Simplex<2> --> Hypercube<2>
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertIndex<Shape::Hypercube<2>, 2>
      {
        static constexpr int shape_dim = 2;
        static constexpr int cell_dim = 2;
        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::Simplex<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioe = index_offsets[1];
          const Index iot = index_offsets[2];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<3> IndexSetTypeIn_v_t;
          typedef IndexSet<3> IndexSetTypeIn_e_t;

          // typedef index vector type
          typedef IndexSetTypeIn_v_t::ConstIndexVectorReference ConstIndexVectorReference_v_t;
          typedef IndexSetTypeIn_e_t::ConstIndexVectorReference ConstIndexVectorReference_e_t;

          // fetch the vertices-at-triangle and edge-at-triangle index set
          const IndexSetTypeIn_v_t& index_set_in_v_t = index_set_holder_in.get_index_set<2,0>();
          const IndexSetTypeIn_e_t& index_set_in_e_t = index_set_holder_in.get_index_set<2,1>();

          // fetch number of coarse mesh triangle
          const Index num_trias = index_set_in_e_t.get_num_entities();

         // loop over all coarse mesh triangles
          for(Index i(0); i < num_trias; ++i)
          {
            // fetch coarse mesh vertices-at-triangle and edge-at-triangle index vector
            ConstIndexVectorReference_v_t v_t = index_set_in_v_t[i];
            ConstIndexVectorReference_e_t e_t = index_set_in_e_t[i];

            // fetch fine mesh vertices-at-quads index vectors
            IndexVectorReference q_0 = index_set_out[offset + 3*i + 0];
            IndexVectorReference q_1 = index_set_out[offset + 3*i + 1];
            IndexVectorReference q_2 = index_set_out[offset + 3*i + 2];

            // start with X and end with the corresponding corner
            // follow the convention in between (wiki)
            q_0[0] = iot + i; // X
            q_0[1] = ioe + e_t[1]; // e_1
            q_0[2] = ioe + e_t[2]; // e_2
            q_0[3] = iov + v_t[0]; // v_0

            q_1[0] = iot + i; // X
            q_1[1] = ioe + e_t[2]; // e_2
            q_1[2] = ioe + e_t[0]; // e_0
            q_1[3] = iov + v_t[1]; // v_1

            q_2[0] = iot + i; // X
            q_2[1] = ioe + e_t[0]; // e_0
            q_2[2] = ioe + e_t[1]; // e_1
            q_2[3] = iov + v_t[2]; // v_2

          }

          // return quads count
          return 3*num_trias;
        }
      };

      /**
       * \brief ShapeConvertIndex implementation for Simplex<3> --> Hypercube<2>
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertIndex<Shape::Hypercube<3>, 2>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 2;
        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::Simplex<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index ioe = index_offsets[1];
          const Index iot = index_offsets[2];
          const Index iotet = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<6> IndexSetTypeIn_e_tet;
          typedef IndexSet<4> IndexSetTypeIn_t_tet;

          // typedef index vector type
          typedef IndexSetTypeIn_e_tet::ConstIndexVectorReference ConstIndexVectorReference_e_tet;
          typedef IndexSetTypeIn_t_tet::ConstIndexVectorReference ConstIndexVectorReference_t_tet;

          // fetch the edge-at-tetra and triangle-at-tetra index set
          const IndexSetTypeIn_e_tet& index_set_in_e_tet = index_set_holder_in.get_index_set<3,1>();
          const IndexSetTypeIn_t_tet& index_set_in_t_tet = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh tetras
          const Index num_tetras = index_set_in_e_tet.get_num_entities();

          // loop over all coarse mesh tetras
          for(Index i(0); i < num_tetras; ++i)
          {
            // fetch coarse mesh edge-at-tetra and triangle-at-tetra index vector
            ConstIndexVectorReference_e_tet e_tet = index_set_in_e_tet[i];
            ConstIndexVectorReference_t_tet t_tet = index_set_in_t_tet[i];

            // fetch fine mesh vertices-at-quad index vectors
            IndexVectorReference q_0 = index_set_out[offset + 6*i + 0];
            IndexVectorReference q_1 = index_set_out[offset + 6*i + 1];
            IndexVectorReference q_2 = index_set_out[offset + 6*i + 2];
            IndexVectorReference q_3 = index_set_out[offset + 6*i + 3];
            IndexVectorReference q_4 = index_set_out[offset + 6*i + 4];
            IndexVectorReference q_5 = index_set_out[offset + 6*i + 5];

            // start with X move to the face with lower index
            // and end up with edge vertex index
            q_0[0] = iotet + i; // X
            q_0[1] = iot + t_tet[2]; // t_2
            q_0[2] = iot + t_tet[3]; // t_3
            q_0[3] = ioe + e_tet[0]; // e_0

            q_1[0] = iotet + i; // X
            q_1[1] = iot + t_tet[1]; // t_1
            q_1[2] = iot + t_tet[3]; // t_3
            q_1[3] = ioe + e_tet[1]; // e_1

            q_2[0] = iotet + i; // X
            q_2[1] = iot + t_tet[1]; // t_1
            q_2[2] = iot + t_tet[2]; // t_2
            q_2[3] = ioe + e_tet[2]; // e_2

            q_3[0] = iotet + i; // X
            q_3[1] = iot + t_tet[0]; // t_0
            q_3[2] = iot + t_tet[3]; // t_3
            q_3[3] = ioe + e_tet[3]; // e_3

            q_4[0] = iotet + i; // X
            q_4[1] = iot + t_tet[0]; // t_0
            q_4[2] = iot + t_tet[2]; // t_2
            q_4[3] = ioe + e_tet[4]; // e_4

            q_5[0] = iotet + i; // X
            q_5[1] = iot + t_tet[0]; // t_0
            q_5[2] = iot + t_tet[1]; // t_1
            q_5[3] = ioe + e_tet[5]; // e_5

          }

          // return quad count
          return 6*num_tetras;
        }
      };

      /**
       * \brief ShapeConvertIndex implementation for Simplex<3> --> Hypercube<3>
       *
       * \author Stefan Wahlers, Constantin Christof
       */

      template<>
      struct ShapeConvertIndex<Shape::Hypercube<3>, 3>
      {
        static constexpr int shape_dim = 3;
        static constexpr int cell_dim = 3;
        typedef Shape::Hypercube<shape_dim> ShapeType;
        typedef Shape::Simplex<shape_dim> OtherShapeType;
        typedef Shape::FaceTraits<ShapeType, cell_dim>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static Index refine(
          IndexSetType& index_set_out,
          const Index offset,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // fetch vertex index offsets
          const Index iov = index_offsets[0];
          const Index ioe = index_offsets[1];
          const Index iot = index_offsets[2];
          const Index iotet = index_offsets[3];

          // typedef index set vector reference for output index set
          typedef IndexSetType::IndexVectorReference IndexVectorReference;

          // typedef index set types
          typedef IndexSet<4> IndexSetTypeIn_v_tet;
          typedef IndexSet<6> IndexSetTypeIn_e_tet;
          typedef IndexSet<4> IndexSetTypeIn_t_tet;

          // typedef index vector type
          typedef IndexSetTypeIn_v_tet::ConstIndexVectorReference ConstIndexVectorReference_v_tet;
          typedef IndexSetTypeIn_e_tet::ConstIndexVectorReference ConstIndexVectorReference_e_tet;
          typedef IndexSetTypeIn_t_tet::ConstIndexVectorReference ConstIndexVectorReference_t_tet;

          // fetch the vertices-at-tetra, edge-at-tetra and triangle-at-tetra index set
          const IndexSetTypeIn_v_tet& index_set_in_v_tet = index_set_holder_in.get_index_set<3,0>();
          const IndexSetTypeIn_e_tet& index_set_in_e_tet = index_set_holder_in.get_index_set<3,1>();
          const IndexSetTypeIn_t_tet& index_set_in_t_tet = index_set_holder_in.get_index_set<3,2>();

          // fetch number of coarse mesh tetras
          const Index num_tetras = index_set_in_e_tet.get_num_entities();

          // loop over all coarse mesh tetras
          for(Index i(0); i < num_tetras; ++i)
          {
            // fetch coarse mesh vertices-at-tetra, edge-at-tetra and triangle-at-tetra index vector
            ConstIndexVectorReference_v_tet v_tet = index_set_in_v_tet[i];
            ConstIndexVectorReference_e_tet e_tet = index_set_in_e_tet[i];
            ConstIndexVectorReference_t_tet t_tet = index_set_in_t_tet[i];

            // fetch fine mesh vertices-at-hexa index vectors
            IndexVectorReference c_0 = index_set_out[offset + 4*i + 0];
            IndexVectorReference c_1 = index_set_out[offset + 4*i + 1];
            IndexVectorReference c_2 = index_set_out[offset + 4*i + 2];
            IndexVectorReference c_3 = index_set_out[offset + 4*i + 3];

            // start with X. move the face with lowest index,
            // move to the next lowest face index and proceed according to
            // the conventions (wiki) and end with the corresponding corner
            c_0[0] = iotet + i; // X
            c_0[1] = iot + t_tet[1]; // t_1
            c_0[2] = iot + t_tet[2]; // t_2
            c_0[3] = ioe + e_tet[2]; // e_2
            c_0[4] = iot + t_tet[3]; // t_3
            c_0[5] = ioe + e_tet[1]; // e_1
            c_0[6] = ioe + e_tet[0]; // e_0
            c_0[7] = iov + v_tet[0]; // v_0

            c_1[0] = iotet + i; // X
            c_1[1] = iot + t_tet[0]; // t_0
            c_1[2] = iot + t_tet[2]; // t_2
            c_1[3] = ioe + e_tet[4]; // e_4
            c_1[4] = iot + t_tet[3]; // t_3
            c_1[5] = ioe + e_tet[3]; // e_3
            c_1[6] = ioe + e_tet[0]; // e_0
            c_1[7] = iov + v_tet[1]; // v_1

            c_2[0] = iotet + i; // X
            c_2[1] = iot + t_tet[0]; // t_0
            c_2[2] = iot + t_tet[1]; // t_1
            c_2[3] = ioe + e_tet[5]; // e_5
            c_2[4] = iot + t_tet[3]; // t_3
            c_2[5] = ioe + e_tet[3]; // e_3
            c_2[6] = ioe + e_tet[1]; // e_1
            c_2[7] = iov + v_tet[2]; // v_2

            c_3[0] = iotet + i; // X
            c_3[1] = iot + t_tet[0]; // t_0
            c_3[2] = iot + t_tet[1]; // t_1
            c_3[3] = ioe + e_tet[5]; // e_5
            c_3[4] = iot + t_tet[2]; // t_2
            c_3[5] = ioe + e_tet[4]; // e_4
            c_3[6] = ioe + e_tet[2]; // e_2
            c_3[7] = iov + v_tet[3]; // v_3
          }

          // return hexa count
          return 4*num_tetras;
        }
      };
      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */

      template<
        typename Shape_,
        int cell_dim_,
        // The following "dummy" argument is necessary for partial specialisation;
        // shape_dim_ *must* coincide with Shape_::dimension !!!
        int shape_dim_ = Shape_::dimension>
      struct ShapeConvertIndexShapeWrapper
      {
        static_assert(shape_dim_ == Shape_::dimension, "invalid shape dimension");
        // the case shape_dim_ = cell_dim_ is specialised below
        static_assert(shape_dim_ > cell_dim_, "invalid cell dimension");
        static_assert(cell_dim_ > 0, "invalid cell dimension");

        typedef Shape_ ShapeType;
        typedef typename OtherShape<Shape_>::Type OtherShapeType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static String name()
        {
          return "ShapeConvertIndexShapeWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ","
            + stringify(shape_dim_) + ">";
        }

        static Index refine(
          IndexSetType& index_set_out,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          typedef typename Shape::FaceTraits<ShapeType, shape_dim_ - 1>::ShapeType FacetType;

          // refine facets
          Index offset = ShapeConvertIndexShapeWrapper<FacetType, cell_dim_>
            ::refine(index_set_out, index_offsets, index_set_holder_in);

          // call index refiner and return new offset
          Index num_faces = ShapeConvertIndex<Shape_, cell_dim_>
            ::refine(index_set_out, offset, index_offsets, index_set_holder_in);

          // validate number of created indices
          Index num_shapes = index_set_holder_in.template get_index_set<shape_dim_, 0>().get_num_entities();
          XASSERTM(num_faces == (ShapeConvertTraits<OtherShapeType, cell_dim_>::count * num_shapes),
            "ShapeConvertIndex output does not match ShapeConvertTraits prediction");

          // calculate new offset and return
          return offset + num_faces;
        }
      };

      template<
        typename Shape_,
        int cell_dim_>
      struct ShapeConvertIndexShapeWrapper<Shape_, cell_dim_, cell_dim_>
      {
        static_assert(cell_dim_ == Shape_::dimension, "invalid shape dimension");
        static_assert(cell_dim_ > 0, "invalid cell dimension");

        typedef Shape_ ShapeType;
        typedef typename OtherShape<Shape_>::Type OtherShapeType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderType;

        static String name()
        {
          return "ShapeConvertIndexShapeWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ","
            + "," + stringify(cell_dim_) + ">";
        }

        static Index refine(
          IndexSetType& index_set_out,
          const Index index_offsets[],
          const IndexSetHolderType& index_set_holder_in)
        {
          // call index refiner
          Index num_faces = ShapeConvertIndex<Shape_, cell_dim_>
            ::refine(index_set_out, 0, index_offsets, index_set_holder_in);

          // validate number of created indices
          Index num_shapes = index_set_holder_in.template get_index_set<cell_dim_, 0>().get_num_entities();
          XASSERTM(num_faces == (ShapeConvertTraits<OtherShapeType, cell_dim_>::count * num_shapes),
            "ShapeConvertIndex output does not match ShapeConvertTraits prediction");

          // return offset
          return num_faces;
        }
      };

      /* ********************************************************************************************************* */

      template<
        typename Shape_,
        int cell_dim_ = Shape_::dimension>
      struct ShapeConvertIndexWrapper
      {
        // the case cell_dim_ = 1 is handled by the partial specialisation below
        static_assert(cell_dim_ > 1, "invalid cell dimension");
        static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");

        typedef Shape_ ShapeType;
        typedef typename OtherShape<Shape_>::Type OtherShapeType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderTypeIn;

        static String name()
        {
          return "ShapeConvertIndexWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ">";
        }

        static void refine(
          IndexSetHolderType& index_set_holder_out,
          const Index num_entities[],
          const IndexSetHolderTypeIn& index_set_holder_in)
        {
          // recursive call of IndexRefineWrapper
          ShapeConvertIndexWrapper<Shape_, cell_dim_ - 1>
            ::refine(index_set_holder_out, num_entities, index_set_holder_in);

          // fetch output index set
          IndexSetType& index_set_out = index_set_holder_out.template get_index_set<cell_dim_, 0>();

          // calculate index offsets
          Index index_offsets[Shape_::dimension+1];
          EntityCounter<ShapeConvertTraits,OtherShapeType, 0>::offset(index_offsets, num_entities);

          // call face wrapper
          ShapeConvertIndexShapeWrapper<Shape_, cell_dim_>
            ::refine(index_set_out, index_offsets, index_set_holder_in);
        }
      };

      template<typename Shape_>
      struct ShapeConvertIndexWrapper<Shape_, 1>
      {
        static_assert(1 <= Shape_::dimension, "invalid shape dimension");

        typedef Shape_ ShapeType;
        typedef typename OtherShape<Shape_>::Type OtherShapeType;
        typedef typename Shape::FaceTraits<ShapeType, 1>::ShapeType CellType;
        typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;
        typedef IndexSetHolder<OtherShapeType> IndexSetHolderTypeIn;

        static String name()
        {
          return "ShapeConvertIndexWrapper<" + Shape_::name() + ",1>";
        }

        static void refine(
          IndexSetHolderType& index_set_holder_out,
          const Index num_entities[],
          const IndexSetHolderTypeIn& index_set_holder_in)
        {
          // fetch output index set
          IndexSetType& index_set_out = index_set_holder_out.template get_index_set<1, 0>();

          // calculate index offsets
          Index index_offsets[Shape_::dimension+1];
          EntityCounter<ShapeConvertTraits,OtherShapeType, 0>::offset(index_offsets, num_entities);

          // call face wrapper
          ShapeConvertIndexShapeWrapper<Shape_, 1>::refine(index_set_out, index_offsets, index_set_holder_in);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_SHAPE_CONVERT_INDEX_HPP
