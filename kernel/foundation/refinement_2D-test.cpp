#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/refinement.hpp>
#include<kernel/foundation/export.hpp>
#include<kernel/archs.hpp>

#include<kernel/foundation/mesh_control.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/export_vtk.hpp>

#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class RefinementTest2D:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    RefinementTest2D(const std::string & tag) :
      TaggedTest<Tag_, Index>("RefinementTest2D<" + tag + ">")
    {
    }

    virtual void run() const override
    {
      /*(0,1) (1,1)
       *  *----*
       *  |    |
       *  |    |
       *  *----*
       *(0,0) (1,0)
       */

      //create attributes for vertex coords
      OT_<Attribute<double, OT_>, std::allocator<Attribute<double, OT_> > > attrs;
      attrs.push_back(Attribute<double, OT_>()); //vertex x-coords
      attrs.push_back(Attribute<double, OT_>()); //vertex y-coords

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(1));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(1));

      /*  2    3
       *  *-1--*
       *  2    |
       *  |    3
       *  *--0-*
       *  0    1
       */

      //creating foundation mesh
      Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_>, OT_> m(0);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);

      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);

      m.add_polytope(pl_face);

      m.add_adjacency(pl_vertex, pl_edge, 0, 0);
      m.add_adjacency(pl_vertex, pl_edge, 0, 2);
      m.add_adjacency(pl_vertex, pl_face, 0, 0);

      m.add_adjacency(pl_vertex, pl_edge, 1, 0);
      m.add_adjacency(pl_vertex, pl_edge, 1, 3);
      m.add_adjacency(pl_vertex, pl_face, 1, 0);

      m.add_adjacency(pl_vertex, pl_edge, 2, 1);
      m.add_adjacency(pl_vertex, pl_edge, 2, 2);
      m.add_adjacency(pl_vertex, pl_face, 2, 0);

      m.add_adjacency(pl_vertex, pl_edge, 3, 1);
      m.add_adjacency(pl_vertex, pl_edge, 3, 3);
      m.add_adjacency(pl_vertex, pl_face, 3, 0);

      Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_>, OT_> m_fine(m);

      OT_<std::shared_ptr<HaloBase<Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> > > > halos;
      halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> >(new Halo<1, PLFace, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m_fine)));
      halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> >(new Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m_fine)));
      halos.at(0)->push_back(0);
      halos.at(1)->push_back(3);

      Refinement<Mem::Main,
                 mrt_standard>::execute(m_fine, &halos, attrs);

/*      Index size_set[3];
      MeshControl<dim_2D>::fill_sizes(m_fine, size_set);

      typedef ConformalMesh<Shape::Hypercube<2> > BaseMeshType;
      BaseMeshType basemesh(size_set);
      MeshControl<dim_2D>::fill_adjacencies(m_fine, basemesh);
      MeshControl<dim_2D>::fill_vertex_sets(m_fine, basemesh, attrs.at(0), attrs.at(1));

      ExportVTK<BaseMeshType> vtkexporter(basemesh);
      vtkexporter.write("test_mesh");
*/
/*      MeshExporter<VTK, Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_>, OT_>, OT_<Attribute<double, OT_>, std::allocator<Attribute<double, OT_> > > > exporter(m_fine, attrs);
      exporter.write("ref_test_2D_result.vtk");
*/
      TEST_CHECK_EQUAL(m_fine.get_topologies().at(ipi_face_vertex).size(), 4ul);
      TEST_CHECK_EQUAL(m_fine.get_topologies().at(ipi_edge_vertex).size(), 12ul);
      TEST_CHECK_EQUAL(m_fine.get_topologies().at(ipi_vertex_edge).size(), 9ul);

      IT_ result_vertices_at_edge_0(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 0));
      TEST_CHECK_EQUAL(result_vertices_at_edge_0.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_0.begin(), result_vertices_at_edge_0.end(), 0) != result_vertices_at_edge_0.end());
      TEST_CHECK(std::find(result_vertices_at_edge_0.begin(), result_vertices_at_edge_0.end(), 4) != result_vertices_at_edge_0.end());
      IT_ result_vertices_at_edge_1(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 1));
      TEST_CHECK_EQUAL(result_vertices_at_edge_1.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_1.begin(), result_vertices_at_edge_1.end(), 2) != result_vertices_at_edge_1.end());
      TEST_CHECK(std::find(result_vertices_at_edge_1.begin(), result_vertices_at_edge_1.end(), 5) != result_vertices_at_edge_1.end());
      IT_ result_vertices_at_edge_2(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 2));
      TEST_CHECK_EQUAL(result_vertices_at_edge_2.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_2.begin(), result_vertices_at_edge_2.end(), 0) != result_vertices_at_edge_2.end());
      TEST_CHECK(std::find(result_vertices_at_edge_2.begin(), result_vertices_at_edge_2.end(), 6) != result_vertices_at_edge_2.end());
      IT_ result_vertices_at_edge_3(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 3));
      TEST_CHECK_EQUAL(result_vertices_at_edge_3.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_3.begin(), result_vertices_at_edge_3.end(), 1) != result_vertices_at_edge_3.end());
      TEST_CHECK(std::find(result_vertices_at_edge_3.begin(), result_vertices_at_edge_3.end(), 7) != result_vertices_at_edge_3.end());
      IT_ result_vertices_at_edge_4(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 4));
      TEST_CHECK_EQUAL(result_vertices_at_edge_4.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_4.begin(), result_vertices_at_edge_4.end(), 4) != result_vertices_at_edge_4.end());
      TEST_CHECK(std::find(result_vertices_at_edge_4.begin(), result_vertices_at_edge_4.end(), 1) != result_vertices_at_edge_4.end());
      IT_ result_vertices_at_edge_5(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 5));
      TEST_CHECK_EQUAL(result_vertices_at_edge_5.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_5.begin(), result_vertices_at_edge_5.end(), 5) != result_vertices_at_edge_5.end());
      TEST_CHECK(std::find(result_vertices_at_edge_5.begin(), result_vertices_at_edge_5.end(), 3) != result_vertices_at_edge_5.end());
      IT_ result_vertices_at_edge_6(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 6));
      TEST_CHECK_EQUAL(result_vertices_at_edge_6.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_6.begin(), result_vertices_at_edge_6.end(), 6) != result_vertices_at_edge_6.end());
      TEST_CHECK(std::find(result_vertices_at_edge_6.begin(), result_vertices_at_edge_6.end(), 2) != result_vertices_at_edge_6.end());
      IT_ result_vertices_at_edge_7(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 7));
      TEST_CHECK_EQUAL(result_vertices_at_edge_7.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_7.begin(), result_vertices_at_edge_7.end(), 7) != result_vertices_at_edge_7.end());
      TEST_CHECK(std::find(result_vertices_at_edge_7.begin(), result_vertices_at_edge_7.end(), 3) != result_vertices_at_edge_7.end());
      IT_ result_vertices_at_edge_8(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 8));
      TEST_CHECK_EQUAL(result_vertices_at_edge_8.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_8.begin(), result_vertices_at_edge_8.end(), 4) != result_vertices_at_edge_8.end());
      TEST_CHECK(std::find(result_vertices_at_edge_8.begin(), result_vertices_at_edge_8.end(), 8) != result_vertices_at_edge_8.end());
      IT_ result_vertices_at_edge_9(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 9));
      TEST_CHECK_EQUAL(result_vertices_at_edge_9.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_9.begin(), result_vertices_at_edge_9.end(), 5) != result_vertices_at_edge_9.end());
      TEST_CHECK(std::find(result_vertices_at_edge_9.begin(), result_vertices_at_edge_9.end(), 8) != result_vertices_at_edge_9.end());
      IT_ result_vertices_at_edge_10(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 10));
      TEST_CHECK_EQUAL(result_vertices_at_edge_10.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_10.begin(), result_vertices_at_edge_10.end(), 6) != result_vertices_at_edge_10.end());
      TEST_CHECK(std::find(result_vertices_at_edge_10.begin(), result_vertices_at_edge_10.end(), 8) != result_vertices_at_edge_10.end());
      IT_ result_vertices_at_edge_11(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 11));
      TEST_CHECK_EQUAL(result_vertices_at_edge_11.size(), 2ul);
      TEST_CHECK(std::find(result_vertices_at_edge_11.begin(), result_vertices_at_edge_11.end(), 7) != result_vertices_at_edge_11.end());
      TEST_CHECK(std::find(result_vertices_at_edge_11.begin(), result_vertices_at_edge_11.end(), 8) != result_vertices_at_edge_11.end());

      IT_ result_vertices_at_face_0(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 0));
      TEST_CHECK_EQUAL(result_vertices_at_face_0.size(), 4ul);
      TEST_CHECK(std::find(result_vertices_at_face_0.begin(), result_vertices_at_face_0.end(), 0) != result_vertices_at_face_0.end());
      TEST_CHECK(std::find(result_vertices_at_face_0.begin(), result_vertices_at_face_0.end(), 4) != result_vertices_at_face_0.end());
      TEST_CHECK(std::find(result_vertices_at_face_0.begin(), result_vertices_at_face_0.end(), 8) != result_vertices_at_face_0.end());
      TEST_CHECK(std::find(result_vertices_at_face_0.begin(), result_vertices_at_face_0.end(), 6) != result_vertices_at_face_0.end());
      IT_ result_vertices_at_face_1(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 1));
      TEST_CHECK_EQUAL(result_vertices_at_face_1.size(), 4ul);
      TEST_CHECK(std::find(result_vertices_at_face_1.begin(), result_vertices_at_face_1.end(), 4) != result_vertices_at_face_1.end());
      TEST_CHECK(std::find(result_vertices_at_face_1.begin(), result_vertices_at_face_1.end(), 1) != result_vertices_at_face_1.end());
      TEST_CHECK(std::find(result_vertices_at_face_1.begin(), result_vertices_at_face_1.end(), 7) != result_vertices_at_face_1.end());
      TEST_CHECK(std::find(result_vertices_at_face_1.begin(), result_vertices_at_face_1.end(), 8) != result_vertices_at_face_1.end());
      IT_ result_vertices_at_face_2(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 2));
      TEST_CHECK_EQUAL(result_vertices_at_face_2.size(), 4ul);
      TEST_CHECK(std::find(result_vertices_at_face_2.begin(), result_vertices_at_face_2.end(), 6) != result_vertices_at_face_2.end());
      TEST_CHECK(std::find(result_vertices_at_face_2.begin(), result_vertices_at_face_2.end(), 8) != result_vertices_at_face_2.end());
      TEST_CHECK(std::find(result_vertices_at_face_2.begin(), result_vertices_at_face_2.end(), 5) != result_vertices_at_face_2.end());
      TEST_CHECK(std::find(result_vertices_at_face_2.begin(), result_vertices_at_face_2.end(), 2) != result_vertices_at_face_2.end());
      IT_ result_vertices_at_face_3(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 3));
      TEST_CHECK_EQUAL(result_vertices_at_face_3.size(), 4ul);
      TEST_CHECK(std::find(result_vertices_at_face_3.begin(), result_vertices_at_face_3.end(), 8) != result_vertices_at_face_3.end());
      TEST_CHECK(std::find(result_vertices_at_face_3.begin(), result_vertices_at_face_3.end(), 7) != result_vertices_at_face_3.end());
      TEST_CHECK(std::find(result_vertices_at_face_3.begin(), result_vertices_at_face_3.end(), 3) != result_vertices_at_face_3.end());
      TEST_CHECK(std::find(result_vertices_at_face_3.begin(), result_vertices_at_face_3.end(), 5) != result_vertices_at_face_3.end());

      IT_ result_faces_at_vertex_0(m_fine.get_adjacent_polytopes(pl_vertex, pl_face, 0));
      TEST_CHECK_EQUAL(result_faces_at_vertex_0.size(), 1ul);
      TEST_CHECK(std::find(result_faces_at_vertex_0.begin(), result_faces_at_vertex_0.end(), 0) != result_faces_at_vertex_0.end());
      IT_ result_faces_at_vertex_1(m_fine.get_adjacent_polytopes(pl_vertex, pl_face, 1));
      TEST_CHECK_EQUAL(result_faces_at_vertex_1.size(), 1ul);
      TEST_CHECK(std::find(result_faces_at_vertex_1.begin(), result_faces_at_vertex_1.end(), 1) != result_faces_at_vertex_1.end());
      IT_ result_faces_at_vertex_2(m_fine.get_adjacent_polytopes(pl_vertex, pl_face, 2));
      TEST_CHECK_EQUAL(result_faces_at_vertex_2.size(), 1ul);
      TEST_CHECK(std::find(result_faces_at_vertex_2.begin(), result_faces_at_vertex_2.end(), 2) != result_faces_at_vertex_2.end());
      IT_ result_faces_at_vertex_3(m_fine.get_adjacent_polytopes(pl_vertex, pl_face, 3));
      TEST_CHECK_EQUAL(result_faces_at_vertex_3.size(), 1ul);
      TEST_CHECK(std::find(result_faces_at_vertex_3.begin(), result_faces_at_vertex_3.end(), 3) != result_faces_at_vertex_3.end());
      IT_ result_faces_at_vertex_4(m_fine.get_adjacent_polytopes(pl_vertex, pl_face, 4));
      TEST_CHECK_EQUAL(result_faces_at_vertex_4.size(), 2ul);
      TEST_CHECK(std::find(result_faces_at_vertex_4.begin(), result_faces_at_vertex_4.end(), 0) != result_faces_at_vertex_4.end());
      TEST_CHECK(std::find(result_faces_at_vertex_4.begin(), result_faces_at_vertex_4.end(), 1) != result_faces_at_vertex_4.end());
      IT_ result_faces_at_vertex_5(m_fine.get_adjacent_polytopes(pl_vertex, pl_face, 5));
      TEST_CHECK_EQUAL(result_faces_at_vertex_5.size(), 2ul);
      TEST_CHECK(std::find(result_faces_at_vertex_5.begin(), result_faces_at_vertex_5.end(), 2) != result_faces_at_vertex_5.end());
      TEST_CHECK(std::find(result_faces_at_vertex_5.begin(), result_faces_at_vertex_5.end(), 3) != result_faces_at_vertex_5.end());
      IT_ result_faces_at_vertex_6(m_fine.get_adjacent_polytopes(pl_vertex, pl_face, 6));
      TEST_CHECK_EQUAL(result_faces_at_vertex_6.size(), 2ul);
      TEST_CHECK(std::find(result_faces_at_vertex_6.begin(), result_faces_at_vertex_6.end(), 0) != result_faces_at_vertex_6.end());
      TEST_CHECK(std::find(result_faces_at_vertex_6.begin(), result_faces_at_vertex_6.end(), 2) != result_faces_at_vertex_6.end());
      IT_ result_faces_at_vertex_7(m_fine.get_adjacent_polytopes(pl_vertex, pl_face, 7));
      TEST_CHECK_EQUAL(result_faces_at_vertex_7.size(), 2ul);
      TEST_CHECK(std::find(result_faces_at_vertex_7.begin(), result_faces_at_vertex_7.end(), 1) != result_faces_at_vertex_7.end());
      TEST_CHECK(std::find(result_faces_at_vertex_7.begin(), result_faces_at_vertex_7.end(), 3) != result_faces_at_vertex_7.end());
      IT_ result_faces_at_vertex_8(m_fine.get_adjacent_polytopes(pl_vertex, pl_face, 8));
      TEST_CHECK_EQUAL(result_faces_at_vertex_8.size(), 4ul);
      TEST_CHECK(std::find(result_faces_at_vertex_8.begin(), result_faces_at_vertex_8.end(), 0) != result_faces_at_vertex_8.end());
      TEST_CHECK(std::find(result_faces_at_vertex_8.begin(), result_faces_at_vertex_8.end(), 1) != result_faces_at_vertex_8.end());
      TEST_CHECK(std::find(result_faces_at_vertex_8.begin(), result_faces_at_vertex_8.end(), 2) != result_faces_at_vertex_8.end());
      TEST_CHECK(std::find(result_faces_at_vertex_8.begin(), result_faces_at_vertex_8.end(), 3) != result_faces_at_vertex_8.end());

      TEST_CHECK_EQUAL(halos.at(0)->get_elements().size(), IndexType_(4));
      TEST_CHECK_EQUAL(halos.at(0)->get_elements().at(0), IndexType_(0));
      TEST_CHECK_EQUAL(halos.at(0)->get_elements().at(1), IndexType_(3));
      TEST_CHECK_EQUAL(halos.at(0)->get_elements().at(2), IndexType_(2));
      TEST_CHECK_EQUAL(halos.at(0)->get_elements().at(3), IndexType_(1));

      TEST_CHECK_EQUAL(halos.at(1)->get_elements().size(), IndexType_(2));
      TEST_CHECK_EQUAL(halos.at(1)->get_elements().at(0), IndexType_(3));
      TEST_CHECK_EQUAL(halos.at(1)->get_elements().at(1), IndexType_(7));

      TEST_CHECK_EQUAL(attrs.at(0).get_data().size(), 9ul);
      TEST_CHECK_EQUAL(attrs.at(1).get_data().size(), 9ul);
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(0).get_data().at(0), double(0.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(1).get_data().at(0), double(0.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(0).get_data().at(1), double(1.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(1).get_data().at(1), double(0.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(0).get_data().at(2), double(0.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(1).get_data().at(2), double(1.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(0).get_data().at(3), double(1.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(1).get_data().at(3), double(1.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(0).get_data().at(4), double(0.5), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(1).get_data().at(4), double(0.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(0).get_data().at(5), double(0.5), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(1).get_data().at(5), double(1.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(0).get_data().at(6), double(0.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(1).get_data().at(6), double(0.5), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(0).get_data().at(7), double(1.0), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(1).get_data().at(7), double(0.5), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(0).get_data().at(8), double(0.5), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(attrs.at(1).get_data().at(8), double(0.5), std::numeric_limits<double>::epsilon());
    }
};
RefinementTest2D<Mem::Main, Index, std::vector, std::vector<Index> > ref_test2_cpu_v_v("std::vector, std::vector");
RefinementTest2D<Mem::Main, Index, std::vector, std::deque<Index> > ref_test2_cpu_v_d("std::vector, std::deque");
RefinementTest2D<Mem::Main, Index, std::deque, std::vector<Index> > ref_test2_cpu_d_v("std::deque, std::vector");
RefinementTest2D<Mem::Main, Index, std::deque, std::deque<Index> > ref_test2_cpu_d_d("std::deque, std::deque");
