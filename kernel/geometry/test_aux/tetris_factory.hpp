#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_TETRIS_FACTORY_HPP
#define KERNEL_GEOMETRY_TEST_AUX_TETRIS_FACTORY_HPP 1

#include <kernel/geometry/conformal_mesh.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {
      template<typename Mesh_>
      class TetrisFactory;

      template<int stride_, typename Coord_>
      class TetrisFactory<ConformalMesh<Shape::Quadrilateral, 2, stride_, Coord_> > :
        public Factory<ConformalMesh<FEAT::Shape::Quadrilateral, 2, stride_, Coord_> >
      {
      public:
        /// mesh type
        typedef ConformalMesh<Shape::Quadrilateral, 2, stride_, Coord_> MeshType;
        typedef Factory<MeshType> BaseClass;

        /// vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// index holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

      public:
        TetrisFactory()
        {
        }

        virtual FEAT::Index get_num_entities(int dim)
        {
          switch(dim)
          {
          case 0:
            return 10;
          case 1:
            return 13;
          case 2:
            return 4;
          default:
            return 0;
          }
        }

        virtual void fill_vertex_set(VertexSetType& vertex_set)
        {
          vertex_set[0][0] = Coord_(1.0);
          vertex_set[0][1] = Coord_(0.0);
          vertex_set[1][0] = Coord_(2.0);
          vertex_set[1][1] = Coord_(0.0);
          vertex_set[2][0] = Coord_(3.0);
          vertex_set[2][1] = Coord_(0.0);
          vertex_set[3][0] = Coord_(0.0);
          vertex_set[3][1] = Coord_(1.0);
          vertex_set[4][0] = Coord_(1.0);
          vertex_set[4][1] = Coord_(1.0);
          vertex_set[5][0] = Coord_(2.0);
          vertex_set[5][1] = Coord_(1.0);
          vertex_set[6][0] = Coord_(3.0);
          vertex_set[6][1] = Coord_(1.0);
          vertex_set[7][0] = Coord_(0.0);
          vertex_set[7][1] = Coord_(2.0);
          vertex_set[8][0] = Coord_(1.0);
          vertex_set[8][1] = Coord_(2.0);
          vertex_set[9][0] = Coord_(2.0);
          vertex_set[9][1] = Coord_(2.0);
        }

        virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
        {
          // set vertices-at-edge
          IndexSet<2>& v_e(index_set_holder.template get_index_set<1,0>());
          v_e( 0, 0) = 0;
          v_e( 0, 1) = 1;
          v_e( 1, 0) = 1;
          v_e( 1, 1) = 2;
          v_e( 2, 0) = 0;
          v_e( 2, 1) = 4;
          v_e( 3, 0) = 1;
          v_e( 3, 1) = 5;
          v_e( 4, 0) = 2;
          v_e( 4, 1) = 6;
          v_e( 5, 0) = 3;
          v_e( 5, 1) = 4;
          v_e( 6, 0) = 4;
          v_e( 6, 1) = 5;
          v_e( 7, 0) = 5;
          v_e( 7, 1) = 6;
          v_e( 8, 0) = 3;
          v_e( 8, 1) = 7;
          v_e( 9, 0) = 4;
          v_e( 9, 1) = 8;
          v_e(10, 0) = 5;
          v_e(10, 1) = 9;
          v_e(11, 0) = 7;
          v_e(11, 1) = 8;
          v_e(12, 0) = 8;
          v_e(12, 1) = 9;

          // set vertices-at-quad
          IndexSet<4>& v_q(index_set_holder.template get_index_set<2,0>());
          v_q(0, 0) = 0;
          v_q(0, 1) = 1;
          v_q(0, 2) = 4;
          v_q(0, 3) = 5;
          v_q(1, 0) = 1;
          v_q(1, 1) = 2;
          v_q(1, 2) = 5;
          v_q(1, 3) = 6;
          v_q(2, 0) = 3;
          v_q(2, 1) = 4;
          v_q(2, 2) = 7;
          v_q(2, 3) = 8;
          v_q(3, 0) = 4;
          v_q(3, 1) = 5;
          v_q(3, 2) = 8;
          v_q(3, 3) = 9;

          // set edges-at-quad
          IndexSet<4>& e_q(index_set_holder.template get_index_set<2,1>());
          e_q(0, 0) = 0;
          e_q(0, 1) = 6;
          e_q(0, 2) = 2;
          e_q(0, 3) = 3;
          e_q(1, 0) = 1;
          e_q(1, 1) = 7;
          e_q(1, 2) = 3;
          e_q(1, 3) = 4;
          e_q(2, 0) = 5;
          e_q(2, 1) = 11;
          e_q(2, 2) = 8;
          e_q(2, 3) = 9;
          e_q(3, 0) = 6;
          e_q(3, 1) = 12;
          e_q(3, 2) = 9;
          e_q(3, 3) = 10;
        }
      };
    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_TEST_AUX_TETRIS_FACTORY_HPP
