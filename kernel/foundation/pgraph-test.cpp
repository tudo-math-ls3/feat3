#include<kernel/base_header.hpp>
#include<test_system/test_system.hpp>

#ifdef FEAT_HAVE_PARMETIS
#include<kernel/archs.hpp>
#include<kernel/foundation/pgraph.hpp>
#include<kernel/geometry/index_calculator.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Foundation;
using namespace FEAT::Geometry;


template<typename Tag_= Mem::Main, typename IndexType_ = Index, typename DataType = float>
class PGraphParmetisTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    PGraphParmetisTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("PGraphParmetisTest<" + tag + ">")
    {
    }

    virtual void run() const override
    {
      Util::Communicator comm( Util::Communicator(MPI_COMM_WORLD) );
      PGraphParmetis pg(2, 1, 2, comm);

      //evoke-test creator function for dual graph
      // *--0--*--1--*
      // 0     1     2
      typedef ConformalMesh<Shape::Hypercube<1> > ConfmeshType1D;
      IndexType_* sizes = new IndexType_[2];
      sizes[0] = 3;
      sizes[1] = 2;
      ConfmeshType1D mesh(sizes);
      typename ConfmeshType1D::template IndexSet<1, 0>::Type& target_vertex_at_edge(mesh.template get_index_set<1, 0>());
      target_vertex_at_edge[0][0] = 0;
      target_vertex_at_edge[0][1] = 1;
      target_vertex_at_edge[1][0] = 1;
      target_vertex_at_edge[1][1] = 2;

      PGraphParmetis dual_graph(mesh, 2, Util::Communicator(MPI_COMM_WORLD));
      auto* result_xadj(dual_graph.get_xadj());
      auto* result_adjncy(dual_graph.get_adjncy());

      TEST_CHECK_EQUAL(dual_graph.get_num_vtx(), 2);
      TEST_CHECK_EQUAL(result_xadj[0], 0);
      TEST_CHECK_EQUAL(result_xadj[1], 2);
      TEST_CHECK_EQUAL(result_xadj[2], 4);
      TEST_CHECK_EQUAL(result_adjncy[0], 0);
      TEST_CHECK_EQUAL(result_adjncy[1], 1);
      TEST_CHECK_EQUAL(result_adjncy[2], 0);
      TEST_CHECK_EQUAL(result_adjncy[3], 1);

      //evoke-test creator function for dual graph
      /*  2    3    5
       *  *-1--*--5-*
       *  2  0 |  1 |
       *  |    3    6
       *  *--0-*--4-*
       *  0    1    4
       */
      typedef ConformalMesh<Shape::Hypercube<2> > ConfmeshType2D;
      IndexType_* sizes2 = new IndexType_[3];
      sizes2[0] = 6;
      sizes2[1] = 7;
      sizes2[2] = 2;
      ConfmeshType2D mesh2(sizes2);
      typename ConfmeshType2D::template IndexSet<1, 0>::Type& target_vertex_at_edge2(mesh2.template get_index_set<1, 0>());
      target_vertex_at_edge2[0][0] = 0;
      target_vertex_at_edge2[0][1] = 1;
      target_vertex_at_edge2[1][0] = 2;
      target_vertex_at_edge2[1][1] = 3;
      target_vertex_at_edge2[2][0] = 0;
      target_vertex_at_edge2[2][1] = 2;
      target_vertex_at_edge2[3][0] = 1;
      target_vertex_at_edge2[3][1] = 3;
      target_vertex_at_edge2[4][0] = 1;
      target_vertex_at_edge2[4][1] = 4;
      target_vertex_at_edge2[5][0] = 3;
      target_vertex_at_edge2[5][1] = 5;
      target_vertex_at_edge2[6][0] = 4;
      target_vertex_at_edge2[6][1] = 5;

      typename ConfmeshType2D::template IndexSet<2, 0>::Type& target_vertex_at_face2(mesh2.template get_index_set<2, 0>());
      target_vertex_at_face2[0][0] = 0;
      target_vertex_at_face2[0][1] = 1;
      target_vertex_at_face2[0][2] = 2;
      target_vertex_at_face2[0][3] = 3;
      target_vertex_at_face2[1][0] = 1;
      target_vertex_at_face2[1][1] = 4;
      target_vertex_at_face2[1][2] = 3;
      target_vertex_at_face2[1][3] = 5;

      //face->edge
      typedef typename ConfmeshType2D::ShapeType face_type;
      typedef typename Shape::FaceTraits<face_type, face_type::dimension - 1>::ShapeType edge_type;

      Geometry::IndexTree<edge_type> it(mesh2.get_num_entities(0));
      it.parse(target_vertex_at_edge2);

      typename ConfmeshType2D::template IndexSet<2, 1>::Type& target_edge_at_face2(mesh2.template get_index_set<2, 1>());
      Geometry::IndexCalculator<face_type, face_type::dimension - 1>::compute(it, target_vertex_at_face2, target_edge_at_face2);

      PGraphParmetis dual_graph2(mesh2, 2, Util::Communicator(MPI_COMM_WORLD));
      auto* result_xadj2(dual_graph2.get_xadj());
      auto* result_adjncy2(dual_graph2.get_adjncy());
      TEST_CHECK_EQUAL(dual_graph2.get_num_vtx(), 2);
      TEST_CHECK_EQUAL(result_xadj2[0], 0);
      TEST_CHECK_EQUAL(result_xadj2[1], 2);
      TEST_CHECK_EQUAL(result_xadj2[2], 4);
      TEST_CHECK_EQUAL(result_adjncy2[0], 0);
      TEST_CHECK_EQUAL(result_adjncy2[1], 1);
      TEST_CHECK_EQUAL(result_adjncy2[2], 0);
      TEST_CHECK_EQUAL(result_adjncy2[3], 1);

      //evoke-test creator function for dual graph
      /*     6    7
       *     *----*
       *    /|   /|
       *  2/ | 3/ |
       *  *----*  *
       *  | /4 | /5
       *  |/   |/
       *  *----*
       *  0    1
       */
      typedef ConformalMesh<Shape::Hypercube<3> > ConfmeshType3D;
      IndexType_* sizes3 = new IndexType_[4];
      sizes3[0] = 8;
      sizes3[1] = 12;
      sizes3[2] = 6;
      sizes3[3] = 1;
      ConfmeshType3D mesh3(sizes3);
      typename ConfmeshType3D::template IndexSet<1, 0>::Type& target_vertex_at_edge3(mesh3.template get_index_set<1, 0>());
      target_vertex_at_edge3[0][0] = 0;
      target_vertex_at_edge3[0][1] = 1;
      target_vertex_at_edge3[1][0] = 2;
      target_vertex_at_edge3[1][1] = 3;
      target_vertex_at_edge3[2][0] = 0;
      target_vertex_at_edge3[2][1] = 2;
      target_vertex_at_edge3[3][0] = 1;
      target_vertex_at_edge3[3][1] = 3;
      target_vertex_at_edge3[4][0] = 1;
      target_vertex_at_edge3[4][1] = 5;
      target_vertex_at_edge3[5][0] = 3;
      target_vertex_at_edge3[5][1] = 7;
      target_vertex_at_edge3[6][0] = 0;
      target_vertex_at_edge3[6][1] = 4;
      target_vertex_at_edge3[7][0] = 2;
      target_vertex_at_edge3[7][1] = 6;
      target_vertex_at_edge3[8][0] = 5;
      target_vertex_at_edge3[8][1] = 4;
      target_vertex_at_edge3[9][0] = 7;
      target_vertex_at_edge3[9][1] = 6;
      target_vertex_at_edge3[10][0] = 5;
      target_vertex_at_edge3[10][1] = 7;
      target_vertex_at_edge3[11][0] = 4;
      target_vertex_at_edge3[11][1] = 6;

      typename ConfmeshType3D::template IndexSet<2, 0>::Type& target_vertex_at_face3(mesh3.template get_index_set<2, 0>());
      target_vertex_at_face3[0][0] = 0;
      target_vertex_at_face3[0][1] = 1;
      target_vertex_at_face3[0][2] = 2;
      target_vertex_at_face3[0][3] = 3;
      target_vertex_at_face3[1][0] = 1;
      target_vertex_at_face3[1][1] = 5;
      target_vertex_at_face3[1][2] = 3;
      target_vertex_at_face3[1][3] = 7;
      target_vertex_at_face3[2][0] = 0;
      target_vertex_at_face3[2][1] = 4;
      target_vertex_at_face3[2][2] = 2;
      target_vertex_at_face3[2][3] = 6;
      target_vertex_at_face3[3][0] = 0;
      target_vertex_at_face3[3][1] = 1;
      target_vertex_at_face3[3][2] = 2;
      target_vertex_at_face3[3][3] = 5;
      target_vertex_at_face3[4][0] = 2;
      target_vertex_at_face3[4][1] = 3;
      target_vertex_at_face3[4][2] = 6;
      target_vertex_at_face3[4][3] = 7;
      target_vertex_at_face3[5][0] = 4;
      target_vertex_at_face3[5][1] = 5;
      target_vertex_at_face3[5][2] = 6;
      target_vertex_at_face3[5][3] = 7;

      typename ConfmeshType3D::template IndexSet<3, 0>::Type& target_vertex_at_poly3(mesh3.template get_index_set<3, 0>());
      target_vertex_at_poly3[0][0] = 0;
      target_vertex_at_poly3[0][1] = 1;
      target_vertex_at_poly3[0][2] = 2;
      target_vertex_at_poly3[0][3] = 3;
      target_vertex_at_poly3[0][4] = 4;
      target_vertex_at_poly3[0][5] = 5;
      target_vertex_at_poly3[0][6] = 6;
      target_vertex_at_poly3[0][7] = 7;

      typedef typename ConfmeshType3D::ShapeType poly_type;
      typedef typename Shape::FaceTraits<poly_type, poly_type::dimension - 1>::ShapeType face_type;
      typedef typename Shape::FaceTraits<poly_type, poly_type::dimension - 2>::ShapeType edge_type;

      Geometry::IndexTree<face_type> it_face(mesh3.get_num_entities(0));
      Geometry::IndexTree<edge_type> it_edge(mesh3.get_num_entities(0));

      it_face.parse(target_vertex_at_face3);
      it_edge.parse(target_vertex_at_edge3);


      typename ConfmeshType3D::template IndexSet<2, 1>::Type& target_edge_at_face3(mesh3.template get_index_set<2, 1>());
      typename ConfmeshType3D::template IndexSet<3, 1>::Type& target_edge_at_poly3(mesh3.template get_index_set<3, 1>());
      typename ConfmeshType3D::template IndexSet<3, 2>::Type& target_face_at_poly3(mesh3.template get_index_set<3, 2>());

      Geometry::IndexCalculator<poly_type, poly_type::dimension - 1>::compute(it_face, target_vertex_at_poly3, target_face_at_poly3);
      Geometry::IndexCalculator<face_type, face_type::dimension - 1>::compute(it_edge, target_vertex_at_face3, target_edge_at_face3);
      Geometry::IndexCalculator<poly_type, poly_type::dimension - 2>::compute(it_edge, target_vertex_at_poly3, target_edge_at_poly3);

      PGraphParmetis dual_graph3(mesh3, 1, Util::Communicator(MPI_COMM_WORLD));
      auto* result_xadj3(dual_graph3.get_xadj());
      auto* result_adjncy3(dual_graph3.get_adjncy());
      TEST_CHECK_EQUAL(dual_graph3.get_num_vtx(), 1);
      TEST_CHECK_EQUAL(result_xadj3[0], 0);
      TEST_CHECK_EQUAL(result_xadj3[1], 1);
      TEST_CHECK_EQUAL(result_adjncy3[0], 0);

      delete[] sizes;
      delete[] sizes2;
      delete[] sizes3;

    }
};
PGraphParmetisTest<> have_parmetis_test("None, Index, float");
#endif
