#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/refinement.hpp>
#include<kernel/foundation/export.hpp>
#include<kernel/archs.hpp>

#include<kernel/foundation/mesh_control.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/cell_sub_set.hpp>
#include <kernel/geometry/export_vtk.hpp>

#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;

template<typename Tag_, typename IndexType_, typename Algo_, template<typename, typename> class OT_, typename IT_>
class RefinementTest3D:
  public TaggedTest<Tag_, IndexType_, Algo_>
{
  public:
    RefinementTest3D(const std::string & tag) :
      TaggedTest<Tag_, Index, Algo_>("RefinementTest3D<" + tag + ">")
    {
    }

    virtual void run() const
    {
      /*
       * (0,1,1)  (1,1,1)
       *      *----*
       *     /    /|
       *(0,1,0)(1,1,0)
       *   *----*  *(1,0,1)
       *   | /  | /
       *   |/   |/
       *   *----*
       *(0,0,0) (1,0,0)
       */

      //create attributes for vertex coords
      OT_<Attribute<double, OT_>, std::allocator<Attribute<double, OT_> > > attrs;
      attrs.push_back(Attribute<double, OT_>()); //vertex x-coords
      attrs.push_back(Attribute<double, OT_>()); //vertex y-coords
      attrs.push_back(Attribute<double, OT_>()); //vertex z-coords

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(0));
      attrs.at(2).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(0));
      attrs.at(2).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(1));
      attrs.at(2).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(1));
      attrs.at(2).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(0));
      attrs.at(2).get_data().push_back(double(1));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(1));
      attrs.at(2).get_data().push_back(double(1));

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(0));
      attrs.at(2).get_data().push_back(double(1));

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(1));
      attrs.at(2).get_data().push_back(double(1));

      /*  2    3
       *  *-1--*
       *  2    |
       *  |    3
       *  *--0-*
       *  0    1
       */

      //creating foundation mesh
      Foundation::Mesh<Dim3D, Foundation::Topology<IndexType_, OT_, IT_>, OT_> m(0);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);

      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);

      m.add_polytope(pl_face);
      m.add_polytope(pl_face);
      m.add_polytope(pl_face);
      m.add_polytope(pl_face);
      m.add_polytope(pl_face);
      m.add_polytope(pl_face);

      m.add_polytope(pl_polyhedron);

      //vertex 0
      m.add_adjacency(pl_vertex, pl_edge, 0, 0);
      m.add_adjacency(pl_vertex, pl_edge, 0, 2);
      m.add_adjacency(pl_vertex, pl_edge, 0, 10);

      m.add_adjacency(pl_vertex, pl_face, 0, 0);
      m.add_adjacency(pl_vertex, pl_face, 0, 3);
      m.add_adjacency(pl_vertex, pl_face, 0, 4);

      m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

      //vertex 1
      m.add_adjacency(pl_vertex, pl_edge, 1, 0);
      m.add_adjacency(pl_vertex, pl_edge, 1, 3);
      m.add_adjacency(pl_vertex, pl_edge, 1, 4);

      m.add_adjacency(pl_vertex, pl_face, 1, 0);
      m.add_adjacency(pl_vertex, pl_face, 1, 2);
      m.add_adjacency(pl_vertex, pl_face, 1, 4);

      m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

      //vertex 2
      m.add_adjacency(pl_vertex, pl_edge, 2, 1);
      m.add_adjacency(pl_vertex, pl_edge, 2, 2);
      m.add_adjacency(pl_vertex, pl_edge, 2, 11);

      m.add_adjacency(pl_vertex, pl_face, 2, 0);
      m.add_adjacency(pl_vertex, pl_face, 2, 3);
      m.add_adjacency(pl_vertex, pl_face, 2, 5);

      m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

      //vertex 3
      m.add_adjacency(pl_vertex, pl_edge, 3, 1);
      m.add_adjacency(pl_vertex, pl_edge, 3, 3);
      m.add_adjacency(pl_vertex, pl_edge, 3, 5);

      m.add_adjacency(pl_vertex, pl_face, 3, 2);
      m.add_adjacency(pl_vertex, pl_face, 3, 0);
      m.add_adjacency(pl_vertex, pl_face, 3, 5);

      m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

      //vertex 4
      m.add_adjacency(pl_vertex, pl_edge, 4, 4);
      m.add_adjacency(pl_vertex, pl_edge, 4, 6);
      m.add_adjacency(pl_vertex, pl_edge, 4, 7);

      m.add_adjacency(pl_vertex, pl_face, 4, 1);
      m.add_adjacency(pl_vertex, pl_face, 4, 2);
      m.add_adjacency(pl_vertex, pl_face, 4, 4);

      m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

      //vertex 5
      m.add_adjacency(pl_vertex, pl_edge, 5, 5);
      m.add_adjacency(pl_vertex, pl_edge, 5, 6);
      m.add_adjacency(pl_vertex, pl_edge, 5, 8);

      m.add_adjacency(pl_vertex, pl_face, 5, 1);
      m.add_adjacency(pl_vertex, pl_face, 5, 2);
      m.add_adjacency(pl_vertex, pl_face, 5, 5);

      m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

      //vertex 6
      m.add_adjacency(pl_vertex, pl_edge, 6, 7);
      m.add_adjacency(pl_vertex, pl_edge, 6, 9);
      m.add_adjacency(pl_vertex, pl_edge, 6, 10);

      m.add_adjacency(pl_vertex, pl_face, 6, 1);
      m.add_adjacency(pl_vertex, pl_face, 6, 3);
      m.add_adjacency(pl_vertex, pl_face, 6, 4);

      m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

      //vertex 7
      m.add_adjacency(pl_vertex, pl_edge, 7, 8);
      m.add_adjacency(pl_vertex, pl_edge, 7, 9);
      m.add_adjacency(pl_vertex, pl_edge, 7, 11);

      m.add_adjacency(pl_vertex, pl_face, 7, 1);
      m.add_adjacency(pl_vertex, pl_face, 7, 3);
      m.add_adjacency(pl_vertex, pl_face, 7, 5);

      m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

      Foundation::Mesh<Dim3D, Foundation::Topology<IndexType_, OT_, IT_>, OT_> m_fine(m);

      OT_<std::shared_ptr<HaloBase<Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, OT_> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, OT_> > > > halos;
      halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, OT_> >(new Halo<1, PLPolyhedron, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>(m_fine)));
      halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, OT_> >(new Halo<0, PLFace, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>(m_fine)));
      halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, OT_> >(new Halo<0, PLEdge, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>(m_fine)));
      halos.at(0)->push_back(0);
      halos.at(1)->push_back(3);
      halos.at(2)->push_back(3);

      Refinement<Mem::Main,
                 Algo::Generic,
                 mrt_standard>::execute(m_fine, &halos, attrs);

/*      Index size_set[4];
      MeshControl<dim_3D>::fill_sizes(m_fine, size_set);

      typedef ConformalMesh<Shape::Hypercube<3> > BaseMeshType;
            std::cout << "TEST3" << std::endl;
      BaseMeshType basemesh(size_set);

      std::cout << "test" << std::endl;
      MeshControl<dim_3D>::fill_adjacencies(m_fine, basemesh);
      MeshControl<dim_3D>::fill_vertex_sets(m_fine, basemesh, attrs.at(0), attrs.at(1),attrs.at(2));

      ExportVTK<BaseMeshType> vtkexporter(basemesh);
      vtkexporter.write("test_3mesh.vtk");
*/
/*      MeshExporter<VTK, Mesh<Dim3D, Foundation::Topology<IndexType_, OT_, IT_>, OT_>, OT_<Attribute<double, OT_>, std::allocator<Attribute<double, OT_> > > > exporter(m_fine, attrs);
      exporter.write("ref_test_3D_result.vtk");
*/
      // right refinement
      IT_ result_vertices_at_edge_0(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 0));

      //flag vector to kick var-tracking limits
      std::vector<bool> res;
      res.push_back(result_vertices_at_edge_0.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_0.begin(), result_vertices_at_edge_0.end(), 0) != result_vertices_at_edge_0.end());
      res.push_back(std::find(result_vertices_at_edge_0.begin(), result_vertices_at_edge_0.end(), 8) != result_vertices_at_edge_0.end());
      IT_ result_vertices_at_edge_1(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 1));
      res.push_back(result_vertices_at_edge_1.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_1.begin(), result_vertices_at_edge_1.end(), 2) != result_vertices_at_edge_1.end());
      res.push_back(std::find(result_vertices_at_edge_1.begin(), result_vertices_at_edge_1.end(), 9) != result_vertices_at_edge_1.end());
      IT_ result_vertices_at_edge_2(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 2));
      res.push_back(result_vertices_at_edge_2.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_2.begin(), result_vertices_at_edge_2.end(), 0) != result_vertices_at_edge_2.end());
      res.push_back(std::find(result_vertices_at_edge_2.begin(), result_vertices_at_edge_2.end(), 10) != result_vertices_at_edge_2.end());
      IT_ result_vertices_at_edge_3(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 3));
      res.push_back(result_vertices_at_edge_3.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_3.begin(), result_vertices_at_edge_3.end(), 1) != result_vertices_at_edge_3.end());
      res.push_back(std::find(result_vertices_at_edge_3.begin(), result_vertices_at_edge_3.end(), 11) != result_vertices_at_edge_3.end());
      IT_ result_vertices_at_edge_4(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 4));
      res.push_back(result_vertices_at_edge_4.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_4.begin(), result_vertices_at_edge_4.end(), 1) != result_vertices_at_edge_4.end());
      res.push_back(std::find(result_vertices_at_edge_4.begin(), result_vertices_at_edge_4.end(), 12) != result_vertices_at_edge_4.end());
      IT_ result_vertices_at_edge_5(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 5));
      res.push_back(result_vertices_at_edge_5.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_5.begin(), result_vertices_at_edge_5.end(), 3) != result_vertices_at_edge_5.end());
      res.push_back(std::find(result_vertices_at_edge_5.begin(), result_vertices_at_edge_5.end(), 13) != result_vertices_at_edge_5.end());
      IT_ result_vertices_at_edge_6(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 6));
      res.push_back(result_vertices_at_edge_6.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_6.begin(), result_vertices_at_edge_6.end(), 14) != result_vertices_at_edge_6.end());
      res.push_back(std::find(result_vertices_at_edge_6.begin(), result_vertices_at_edge_6.end(), 4) != result_vertices_at_edge_6.end());
      IT_ result_vertices_at_edge_7(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 7));
      res.push_back(result_vertices_at_edge_7.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_7.begin(), result_vertices_at_edge_7.end(), 15) != result_vertices_at_edge_7.end());
      res.push_back(std::find(result_vertices_at_edge_7.begin(), result_vertices_at_edge_7.end(), 4) != result_vertices_at_edge_7.end());
      IT_ result_vertices_at_edge_8(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 8));
      res.push_back(result_vertices_at_edge_8.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_8.begin(), result_vertices_at_edge_8.end(), 16) != result_vertices_at_edge_8.end());
      res.push_back(std::find(result_vertices_at_edge_8.begin(), result_vertices_at_edge_8.end(), 5) != result_vertices_at_edge_8.end());
      IT_ result_vertices_at_edge_9(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 9));
      res.push_back(result_vertices_at_edge_9.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_9.begin(), result_vertices_at_edge_9.end(), 17) != result_vertices_at_edge_9.end());
      res.push_back(std::find(result_vertices_at_edge_9.begin(), result_vertices_at_edge_9.end(), 6) != result_vertices_at_edge_9.end());
      IT_ result_vertices_at_edge_10(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 10));
      res.push_back(result_vertices_at_edge_10.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_10.begin(), result_vertices_at_edge_10.end(), 0) != result_vertices_at_edge_10.end());
      res.push_back(std::find(result_vertices_at_edge_10.begin(), result_vertices_at_edge_10.end(), 18) != result_vertices_at_edge_10.end());
      IT_ result_vertices_at_edge_11(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 11));
      res.push_back(result_vertices_at_edge_11.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_11.begin(), result_vertices_at_edge_11.end(), 2) != result_vertices_at_edge_11.end());
      res.push_back(std::find(result_vertices_at_edge_11.begin(), result_vertices_at_edge_11.end(), 19) != result_vertices_at_edge_11.end());
      IT_ result_vertices_at_edge_12(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 12));
      res.push_back(result_vertices_at_edge_12.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_12.begin(), result_vertices_at_edge_12.end(), 8) != result_vertices_at_edge_12.end());
      res.push_back(std::find(result_vertices_at_edge_12.begin(), result_vertices_at_edge_12.end(), 1) != result_vertices_at_edge_12.end());
      IT_ result_vertices_at_edge_13(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 13));
      res.push_back(result_vertices_at_edge_13.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_13.begin(), result_vertices_at_edge_13.end(), 9) != result_vertices_at_edge_13.end());
      res.push_back(std::find(result_vertices_at_edge_13.begin(), result_vertices_at_edge_13.end(), 3) != result_vertices_at_edge_13.end());
      IT_ result_vertices_at_edge_14(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 14));
      res.push_back(result_vertices_at_edge_14.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_14.begin(), result_vertices_at_edge_14.end(), 2) != result_vertices_at_edge_14.end());
      res.push_back(std::find(result_vertices_at_edge_14.begin(), result_vertices_at_edge_14.end(), 10) != result_vertices_at_edge_14.end());
      IT_ result_vertices_at_edge_15(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 15));
      res.push_back(result_vertices_at_edge_15.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_15.begin(), result_vertices_at_edge_15.end(), 3) != result_vertices_at_edge_15.end());
      res.push_back(std::find(result_vertices_at_edge_15.begin(), result_vertices_at_edge_15.end(), 11) != result_vertices_at_edge_15.end());
      IT_ result_vertices_at_edge_16(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 16));
      res.push_back(result_vertices_at_edge_16.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_16.begin(), result_vertices_at_edge_16.end(), 12) != result_vertices_at_edge_16.end());
      res.push_back(std::find(result_vertices_at_edge_16.begin(), result_vertices_at_edge_16.end(), 4) != result_vertices_at_edge_16.end());
      IT_ result_vertices_at_edge_17(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 17));
      res.push_back(result_vertices_at_edge_17.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_17.begin(), result_vertices_at_edge_17.end(), 5) != result_vertices_at_edge_17.end());
      res.push_back(std::find(result_vertices_at_edge_17.begin(), result_vertices_at_edge_17.end(), 13) != result_vertices_at_edge_17.end());
      IT_ result_vertices_at_edge_18(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 18));
      res.push_back(result_vertices_at_edge_18.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_18.begin(), result_vertices_at_edge_18.end(), 5) != result_vertices_at_edge_18.end());
      res.push_back(std::find(result_vertices_at_edge_18.begin(), result_vertices_at_edge_18.end(), 14) != result_vertices_at_edge_18.end());
      IT_ result_vertices_at_edge_19(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 19));
      res.push_back(result_vertices_at_edge_19.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_19.begin(), result_vertices_at_edge_19.end(), 6) != result_vertices_at_edge_19.end());
      res.push_back(std::find(result_vertices_at_edge_19.begin(), result_vertices_at_edge_19.end(), 15) != result_vertices_at_edge_19.end());
      IT_ result_vertices_at_edge_20(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 20));
      res.push_back(result_vertices_at_edge_20.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_20.begin(), result_vertices_at_edge_20.end(), 7) != result_vertices_at_edge_20.end());
      res.push_back(std::find(result_vertices_at_edge_20.begin(), result_vertices_at_edge_20.end(), 16) != result_vertices_at_edge_20.end());
      IT_ result_vertices_at_edge_21(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 21));
      res.push_back(result_vertices_at_edge_21.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_21.begin(), result_vertices_at_edge_21.end(), 7) != result_vertices_at_edge_21.end());
      res.push_back(std::find(result_vertices_at_edge_21.begin(), result_vertices_at_edge_21.end(), 17) != result_vertices_at_edge_21.end());
      IT_ result_vertices_at_edge_22(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 22));
      res.push_back(result_vertices_at_edge_22.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_22.begin(), result_vertices_at_edge_22.end(), 6) != result_vertices_at_edge_22.end());
      res.push_back(std::find(result_vertices_at_edge_22.begin(), result_vertices_at_edge_22.end(), 18) != result_vertices_at_edge_22.end());
      IT_ result_vertices_at_edge_23(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 23));
      res.push_back(result_vertices_at_edge_23.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_23.begin(), result_vertices_at_edge_23.end(), 7) != result_vertices_at_edge_23.end());
      res.push_back(std::find(result_vertices_at_edge_23.begin(), result_vertices_at_edge_23.end(), 19) != result_vertices_at_edge_23.end());
      IT_ result_vertices_at_edge_24(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 24));
      res.push_back(result_vertices_at_edge_24.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_24.begin(), result_vertices_at_edge_24.end(), 20) != result_vertices_at_edge_24.end());
      res.push_back(std::find(result_vertices_at_edge_24.begin(), result_vertices_at_edge_24.end(), 8) != result_vertices_at_edge_24.end());
      IT_ result_vertices_at_edge_25(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 25));
      res.push_back(result_vertices_at_edge_25.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_25.begin(), result_vertices_at_edge_25.end(), 9) != result_vertices_at_edge_25.end());
      res.push_back(std::find(result_vertices_at_edge_25.begin(), result_vertices_at_edge_25.end(), 20) != result_vertices_at_edge_25.end());
      IT_ result_vertices_at_edge_26(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 26));
      res.push_back(result_vertices_at_edge_26.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_26.begin(), result_vertices_at_edge_26.end(), 10) != result_vertices_at_edge_26.end());
      res.push_back(std::find(result_vertices_at_edge_26.begin(), result_vertices_at_edge_26.end(), 20) != result_vertices_at_edge_26.end());
      IT_ result_vertices_at_edge_27(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 27));
      res.push_back(result_vertices_at_edge_27.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_27.begin(), result_vertices_at_edge_27.end(), 20) != result_vertices_at_edge_27.end());
      res.push_back(std::find(result_vertices_at_edge_27.begin(), result_vertices_at_edge_27.end(), 11) != result_vertices_at_edge_27.end());
      IT_ result_vertices_at_edge_28(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 28));
      res.push_back(result_vertices_at_edge_28.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_28.begin(), result_vertices_at_edge_28.end(), 21) != result_vertices_at_edge_28.end());
      res.push_back(std::find(result_vertices_at_edge_28.begin(), result_vertices_at_edge_28.end(), 14) != result_vertices_at_edge_28.end());
      IT_ result_vertices_at_edge_29(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 29));
      res.push_back(result_vertices_at_edge_29.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_29.begin(), result_vertices_at_edge_29.end(), 21) != result_vertices_at_edge_29.end());
      res.push_back(std::find(result_vertices_at_edge_29.begin(), result_vertices_at_edge_29.end(), 15) != result_vertices_at_edge_29.end());
      IT_ result_vertices_at_edge_30(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 30));
      res.push_back(result_vertices_at_edge_30.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_30.begin(), result_vertices_at_edge_30.end(), 16) != result_vertices_at_edge_30.end());
      res.push_back(std::find(result_vertices_at_edge_30.begin(), result_vertices_at_edge_30.end(), 21) != result_vertices_at_edge_30.end());
      IT_ result_vertices_at_edge_31(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 31));
      res.push_back(result_vertices_at_edge_31.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_31.begin(), result_vertices_at_edge_31.end(), 17) != result_vertices_at_edge_31.end());
      res.push_back(std::find(result_vertices_at_edge_31.begin(), result_vertices_at_edge_31.end(), 21) != result_vertices_at_edge_31.end());
      IT_ result_vertices_at_edge_32(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 32));
      res.push_back(result_vertices_at_edge_32.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_32.begin(), result_vertices_at_edge_32.end(), 11) != result_vertices_at_edge_32.end());
      res.push_back(std::find(result_vertices_at_edge_32.begin(), result_vertices_at_edge_32.end(), 22) != result_vertices_at_edge_32.end());
      IT_ result_vertices_at_edge_33(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 33));
      res.push_back(result_vertices_at_edge_33.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_33.begin(), result_vertices_at_edge_33.end(), 22) != result_vertices_at_edge_33.end());
      res.push_back(std::find(result_vertices_at_edge_33.begin(), result_vertices_at_edge_33.end(), 12) != result_vertices_at_edge_33.end());
      IT_ result_vertices_at_edge_34(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 34));
      res.push_back(result_vertices_at_edge_34.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_34.begin(), result_vertices_at_edge_34.end(), 13) != result_vertices_at_edge_34.end());
      res.push_back(std::find(result_vertices_at_edge_34.begin(), result_vertices_at_edge_34.end(), 22) != result_vertices_at_edge_34.end());
      IT_ result_vertices_at_edge_35(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 35));
      res.push_back(result_vertices_at_edge_35.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_35.begin(), result_vertices_at_edge_35.end(), 22) != result_vertices_at_edge_35.end());
      res.push_back(std::find(result_vertices_at_edge_35.begin(), result_vertices_at_edge_35.end(), 14) != result_vertices_at_edge_35.end());
      IT_ result_vertices_at_edge_36(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 36));
      res.push_back(result_vertices_at_edge_36.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_36.begin(), result_vertices_at_edge_36.end(), 10) != result_vertices_at_edge_36.end());
      res.push_back(std::find(result_vertices_at_edge_36.begin(), result_vertices_at_edge_36.end(), 23) != result_vertices_at_edge_36.end());
      IT_ result_vertices_at_edge_37(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 37));
      res.push_back(result_vertices_at_edge_37.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_37.begin(), result_vertices_at_edge_37.end(), 17) != result_vertices_at_edge_37.end());
      res.push_back(std::find(result_vertices_at_edge_37.begin(), result_vertices_at_edge_37.end(), 23) != result_vertices_at_edge_37.end());
      IT_ result_vertices_at_edge_38(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 38));
      res.push_back(result_vertices_at_edge_38.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_38.begin(), result_vertices_at_edge_38.end(), 23) != result_vertices_at_edge_38.end());
      res.push_back(std::find(result_vertices_at_edge_38.begin(), result_vertices_at_edge_38.end(), 18) != result_vertices_at_edge_38.end());
      IT_ result_vertices_at_edge_39(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 39));
      res.push_back(result_vertices_at_edge_39.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_39.begin(), result_vertices_at_edge_39.end(), 19) != result_vertices_at_edge_39.end());
      res.push_back(std::find(result_vertices_at_edge_39.begin(), result_vertices_at_edge_39.end(), 23) != result_vertices_at_edge_39.end());
      IT_ result_vertices_at_edge_40(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 40));
      res.push_back(result_vertices_at_edge_40.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_40.begin(), result_vertices_at_edge_40.end(), 8) != result_vertices_at_edge_40.end());
      res.push_back(std::find(result_vertices_at_edge_40.begin(), result_vertices_at_edge_40.end(), 24) != result_vertices_at_edge_40.end());
      IT_ result_vertices_at_edge_41(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 41));
      res.push_back(result_vertices_at_edge_41.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_41.begin(), result_vertices_at_edge_41.end(), 24) != result_vertices_at_edge_41.end());
      res.push_back(std::find(result_vertices_at_edge_41.begin(), result_vertices_at_edge_41.end(), 12) != result_vertices_at_edge_41.end());
      IT_ result_vertices_at_edge_42(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 42));
      res.push_back(result_vertices_at_edge_42.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_42.begin(), result_vertices_at_edge_42.end(), 24) != result_vertices_at_edge_42.end());
      res.push_back(std::find(result_vertices_at_edge_42.begin(), result_vertices_at_edge_42.end(), 15) != result_vertices_at_edge_42.end());
      IT_ result_vertices_at_edge_43(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 43));
      res.push_back(result_vertices_at_edge_43.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_43.begin(), result_vertices_at_edge_43.end(), 18) != result_vertices_at_edge_43.end());
      res.push_back(std::find(result_vertices_at_edge_43.begin(), result_vertices_at_edge_43.end(), 24) != result_vertices_at_edge_43.end());
      IT_ result_vertices_at_edge_44(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 44));
      res.push_back(result_vertices_at_edge_44.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_44.begin(), result_vertices_at_edge_44.end(), 9) != result_vertices_at_edge_44.end());
      res.push_back(std::find(result_vertices_at_edge_44.begin(), result_vertices_at_edge_44.end(), 25) != result_vertices_at_edge_44.end());
      IT_ result_vertices_at_edge_45(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 45));
      res.push_back(result_vertices_at_edge_45.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_45.begin(), result_vertices_at_edge_45.end(), 25) != result_vertices_at_edge_45.end());
      res.push_back(std::find(result_vertices_at_edge_45.begin(), result_vertices_at_edge_45.end(), 13) != result_vertices_at_edge_45.end());
      IT_ result_vertices_at_edge_46(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 46));
      res.push_back(result_vertices_at_edge_46.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_46.begin(), result_vertices_at_edge_46.end(), 25) != result_vertices_at_edge_46.end());
      res.push_back(std::find(result_vertices_at_edge_46.begin(), result_vertices_at_edge_46.end(), 16) != result_vertices_at_edge_46.end());
      IT_ result_vertices_at_edge_47(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 47));
      res.push_back(result_vertices_at_edge_47.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_47.begin(), result_vertices_at_edge_47.end(), 19) != result_vertices_at_edge_47.end());
      res.push_back(std::find(result_vertices_at_edge_47.begin(), result_vertices_at_edge_47.end(), 25) != result_vertices_at_edge_47.end());
      IT_ result_vertices_at_edge_48(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 48));
      res.push_back(result_vertices_at_edge_48.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_48.begin(), result_vertices_at_edge_48.end(), 20) != result_vertices_at_edge_48.end());
      res.push_back(std::find(result_vertices_at_edge_48.begin(), result_vertices_at_edge_48.end(), 26) != result_vertices_at_edge_48.end());
      IT_ result_vertices_at_edge_49(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 49));
      res.push_back(result_vertices_at_edge_49.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_49.begin(), result_vertices_at_edge_49.end(), 26) != result_vertices_at_edge_49.end());
      res.push_back(std::find(result_vertices_at_edge_49.begin(), result_vertices_at_edge_49.end(), 21) != result_vertices_at_edge_49.end());
      IT_ result_vertices_at_edge_50(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 50));
      res.push_back(result_vertices_at_edge_50.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_50.begin(), result_vertices_at_edge_50.end(), 26) != result_vertices_at_edge_50.end());
      res.push_back(std::find(result_vertices_at_edge_50.begin(), result_vertices_at_edge_50.end(), 22) != result_vertices_at_edge_50.end());
      IT_ result_vertices_at_edge_51(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 51));
      res.push_back(result_vertices_at_edge_51.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_51.begin(), result_vertices_at_edge_51.end(), 23) != result_vertices_at_edge_51.end());
      res.push_back(std::find(result_vertices_at_edge_51.begin(), result_vertices_at_edge_51.end(), 26) != result_vertices_at_edge_51.end());
      IT_ result_vertices_at_edge_52(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 52));
      res.push_back(result_vertices_at_edge_52.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_52.begin(), result_vertices_at_edge_52.end(), 24) != result_vertices_at_edge_52.end());
      res.push_back(std::find(result_vertices_at_edge_52.begin(), result_vertices_at_edge_52.end(), 26) != result_vertices_at_edge_52.end());
      IT_ result_vertices_at_edge_53(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 53));
      res.push_back(result_vertices_at_edge_53.size() ==  2ul);
      res.push_back(std::find(result_vertices_at_edge_53.begin(), result_vertices_at_edge_53.end(), 25) != result_vertices_at_edge_53.end());
      res.push_back(std::find(result_vertices_at_edge_53.begin(), result_vertices_at_edge_53.end(), 26) != result_vertices_at_edge_53.end());

      IT_ result_vertices_at_face_0(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 0));
      res.push_back(result_vertices_at_face_0.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_0.begin(), result_vertices_at_face_0.end(), 0) != result_vertices_at_face_0.end());
      res.push_back(std::find(result_vertices_at_face_0.begin(), result_vertices_at_face_0.end(), 8) != result_vertices_at_face_0.end());
      res.push_back(std::find(result_vertices_at_face_0.begin(), result_vertices_at_face_0.end(), 10) != result_vertices_at_face_0.end());
      res.push_back(std::find(result_vertices_at_face_0.begin(), result_vertices_at_face_0.end(), 20) != result_vertices_at_face_0.end());
      IT_ result_vertices_at_face_1(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 1));
      res.push_back(result_vertices_at_face_1.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_1.begin(), result_vertices_at_face_1.end(), 15) != result_vertices_at_face_1.end());
      res.push_back(std::find(result_vertices_at_face_1.begin(), result_vertices_at_face_1.end(), 4) != result_vertices_at_face_1.end());
      res.push_back(std::find(result_vertices_at_face_1.begin(), result_vertices_at_face_1.end(), 14) != result_vertices_at_face_1.end());
      res.push_back(std::find(result_vertices_at_face_1.begin(), result_vertices_at_face_1.end(), 21) != result_vertices_at_face_1.end());
      IT_ result_vertices_at_face_2(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 2));
      res.push_back(result_vertices_at_face_2.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_2.begin(), result_vertices_at_face_2.end(), 1) != result_vertices_at_face_2.end());
      res.push_back(std::find(result_vertices_at_face_2.begin(), result_vertices_at_face_2.end(), 12) != result_vertices_at_face_2.end());
      res.push_back(std::find(result_vertices_at_face_2.begin(), result_vertices_at_face_2.end(), 22) != result_vertices_at_face_2.end());
      res.push_back(std::find(result_vertices_at_face_2.begin(), result_vertices_at_face_2.end(), 11) != result_vertices_at_face_2.end());
      IT_ result_vertices_at_face_3(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 3));
      res.push_back(result_vertices_at_face_3.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_3.begin(), result_vertices_at_face_3.end(), 0) != result_vertices_at_face_3.end());
      res.push_back(std::find(result_vertices_at_face_3.begin(), result_vertices_at_face_3.end(), 18) != result_vertices_at_face_3.end());
      res.push_back(std::find(result_vertices_at_face_3.begin(), result_vertices_at_face_3.end(), 23) != result_vertices_at_face_3.end());
      res.push_back(std::find(result_vertices_at_face_3.begin(), result_vertices_at_face_3.end(), 10) != result_vertices_at_face_3.end());
      IT_ result_vertices_at_face_4(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 4));
      res.push_back(result_vertices_at_face_4.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_4.begin(), result_vertices_at_face_4.end(), 0) != result_vertices_at_face_4.end());
      res.push_back(std::find(result_vertices_at_face_4.begin(), result_vertices_at_face_4.end(), 8) != result_vertices_at_face_4.end());
      res.push_back(std::find(result_vertices_at_face_4.begin(), result_vertices_at_face_4.end(), 24) != result_vertices_at_face_4.end());
      res.push_back(std::find(result_vertices_at_face_4.begin(), result_vertices_at_face_4.end(), 18) != result_vertices_at_face_4.end());
      IT_ result_vertices_at_face_5(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 5));
      res.push_back(result_vertices_at_face_5.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_5.begin(), result_vertices_at_face_5.end(), 2) != result_vertices_at_face_5.end());
      res.push_back(std::find(result_vertices_at_face_5.begin(), result_vertices_at_face_5.end(), 9) != result_vertices_at_face_5.end());
      res.push_back(std::find(result_vertices_at_face_5.begin(), result_vertices_at_face_5.end(), 25) != result_vertices_at_face_5.end());
      res.push_back(std::find(result_vertices_at_face_5.begin(), result_vertices_at_face_5.end(), 19) != result_vertices_at_face_5.end());
      IT_ result_vertices_at_face_6(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 6));
      res.push_back(result_vertices_at_face_6.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_6.begin(), result_vertices_at_face_6.end(), 8) != result_vertices_at_face_6.end());
      res.push_back(std::find(result_vertices_at_face_6.begin(), result_vertices_at_face_6.end(), 1) != result_vertices_at_face_6.end());
      res.push_back(std::find(result_vertices_at_face_6.begin(), result_vertices_at_face_6.end(), 11) != result_vertices_at_face_6.end());
      res.push_back(std::find(result_vertices_at_face_6.begin(), result_vertices_at_face_6.end(), 20) != result_vertices_at_face_6.end());
      IT_ result_vertices_at_face_7(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 7));
      res.push_back(result_vertices_at_face_7.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_7.begin(), result_vertices_at_face_7.end(), 10) != result_vertices_at_face_7.end());
      res.push_back(std::find(result_vertices_at_face_7.begin(), result_vertices_at_face_7.end(), 20) != result_vertices_at_face_7.end());
      res.push_back(std::find(result_vertices_at_face_7.begin(), result_vertices_at_face_7.end(), 9) != result_vertices_at_face_7.end());
      res.push_back(std::find(result_vertices_at_face_7.begin(), result_vertices_at_face_7.end(), 2) != result_vertices_at_face_7.end());
      IT_ result_vertices_at_face_8(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 8));
      res.push_back(result_vertices_at_face_8.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_8.begin(), result_vertices_at_face_8.end(), 20) != result_vertices_at_face_8.end());
      res.push_back(std::find(result_vertices_at_face_8.begin(), result_vertices_at_face_8.end(), 11) != result_vertices_at_face_8.end());
      res.push_back(std::find(result_vertices_at_face_8.begin(), result_vertices_at_face_8.end(), 3) != result_vertices_at_face_8.end());
      res.push_back(std::find(result_vertices_at_face_8.begin(), result_vertices_at_face_8.end(), 9) != result_vertices_at_face_8.end());
      IT_ result_vertices_at_face_9(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 9));
      res.push_back(result_vertices_at_face_9.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_9.begin(), result_vertices_at_face_9.end(), 21) != result_vertices_at_face_9.end());
      res.push_back(std::find(result_vertices_at_face_9.begin(), result_vertices_at_face_9.end(), 14) != result_vertices_at_face_9.end());
      res.push_back(std::find(result_vertices_at_face_9.begin(), result_vertices_at_face_9.end(), 5) != result_vertices_at_face_9.end());
      res.push_back(std::find(result_vertices_at_face_9.begin(), result_vertices_at_face_9.end(), 16) != result_vertices_at_face_9.end());
      IT_ result_vertices_at_face_10(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 10));
      res.push_back(result_vertices_at_face_10.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_10.begin(), result_vertices_at_face_10.end(), 6) != result_vertices_at_face_10.end());
      res.push_back(std::find(result_vertices_at_face_10.begin(), result_vertices_at_face_10.end(), 15) != result_vertices_at_face_10.end());
      res.push_back(std::find(result_vertices_at_face_10.begin(), result_vertices_at_face_10.end(), 21) != result_vertices_at_face_10.end());
      res.push_back(std::find(result_vertices_at_face_10.begin(), result_vertices_at_face_10.end(), 17) != result_vertices_at_face_10.end());
      IT_ result_vertices_at_face_11(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 11));
      res.push_back(result_vertices_at_face_11.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_11.begin(), result_vertices_at_face_11.end(), 17) != result_vertices_at_face_11.end());
      res.push_back(std::find(result_vertices_at_face_11.begin(), result_vertices_at_face_11.end(), 21) != result_vertices_at_face_11.end());
      res.push_back(std::find(result_vertices_at_face_11.begin(), result_vertices_at_face_11.end(), 16) != result_vertices_at_face_11.end());
      res.push_back(std::find(result_vertices_at_face_11.begin(), result_vertices_at_face_11.end(), 7) != result_vertices_at_face_11.end());
      IT_ result_vertices_at_face_12(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 12));
      res.push_back(result_vertices_at_face_12.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_12.begin(), result_vertices_at_face_12.end(), 11) != result_vertices_at_face_12.end());
      res.push_back(std::find(result_vertices_at_face_12.begin(), result_vertices_at_face_12.end(), 22) != result_vertices_at_face_12.end());
      res.push_back(std::find(result_vertices_at_face_12.begin(), result_vertices_at_face_12.end(), 13) != result_vertices_at_face_12.end());
      res.push_back(std::find(result_vertices_at_face_12.begin(), result_vertices_at_face_12.end(), 3) != result_vertices_at_face_12.end());
      IT_ result_vertices_at_face_13(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 13));
      res.push_back(result_vertices_at_face_13.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_13.begin(), result_vertices_at_face_13.end(), 12) != result_vertices_at_face_13.end());
      res.push_back(std::find(result_vertices_at_face_13.begin(), result_vertices_at_face_13.end(), 4) != result_vertices_at_face_13.end());
      res.push_back(std::find(result_vertices_at_face_13.begin(), result_vertices_at_face_13.end(), 14) != result_vertices_at_face_13.end());
      res.push_back(std::find(result_vertices_at_face_13.begin(), result_vertices_at_face_13.end(), 22) != result_vertices_at_face_13.end());
      IT_ result_vertices_at_face_14(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 14));
      res.push_back(result_vertices_at_face_14.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_14.begin(), result_vertices_at_face_14.end(), 22) != result_vertices_at_face_14.end());
      res.push_back(std::find(result_vertices_at_face_14.begin(), result_vertices_at_face_14.end(), 14) != result_vertices_at_face_14.end());
      res.push_back(std::find(result_vertices_at_face_14.begin(), result_vertices_at_face_14.end(), 5) != result_vertices_at_face_14.end());
      res.push_back(std::find(result_vertices_at_face_14.begin(), result_vertices_at_face_14.end(), 13) != result_vertices_at_face_14.end());
      IT_ result_vertices_at_face_15(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 15));
      res.push_back(result_vertices_at_face_15.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_15.begin(), result_vertices_at_face_15.end(), 10) != result_vertices_at_face_15.end());
      res.push_back(std::find(result_vertices_at_face_15.begin(), result_vertices_at_face_15.end(), 23) != result_vertices_at_face_15.end());
      res.push_back(std::find(result_vertices_at_face_15.begin(), result_vertices_at_face_15.end(), 19) != result_vertices_at_face_15.end());
      res.push_back(std::find(result_vertices_at_face_15.begin(), result_vertices_at_face_15.end(), 2) != result_vertices_at_face_15.end());
      IT_ result_vertices_at_face_16(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 16));
      res.push_back(result_vertices_at_face_16.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_16.begin(), result_vertices_at_face_16.end(), 18) != result_vertices_at_face_16.end());
      res.push_back(std::find(result_vertices_at_face_16.begin(), result_vertices_at_face_16.end(), 6) != result_vertices_at_face_16.end());
      res.push_back(std::find(result_vertices_at_face_16.begin(), result_vertices_at_face_16.end(), 17) != result_vertices_at_face_16.end());
      res.push_back(std::find(result_vertices_at_face_16.begin(), result_vertices_at_face_16.end(), 23) != result_vertices_at_face_16.end());
      IT_ result_vertices_at_face_17(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 17));
      res.push_back(result_vertices_at_face_17.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_17.begin(), result_vertices_at_face_17.end(), 23) != result_vertices_at_face_17.end());
      res.push_back(std::find(result_vertices_at_face_17.begin(), result_vertices_at_face_17.end(), 17) != result_vertices_at_face_17.end());
      res.push_back(std::find(result_vertices_at_face_17.begin(), result_vertices_at_face_17.end(), 7) != result_vertices_at_face_17.end());
      res.push_back(std::find(result_vertices_at_face_17.begin(), result_vertices_at_face_17.end(), 19) != result_vertices_at_face_17.end());
      IT_ result_vertices_at_face_18(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 18));
      res.push_back(result_vertices_at_face_18.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_18.begin(), result_vertices_at_face_18.end(), 8) != result_vertices_at_face_18.end());
      res.push_back(std::find(result_vertices_at_face_18.begin(), result_vertices_at_face_18.end(), 1) != result_vertices_at_face_18.end());
      res.push_back(std::find(result_vertices_at_face_18.begin(), result_vertices_at_face_18.end(), 12) != result_vertices_at_face_18.end());
      res.push_back(std::find(result_vertices_at_face_18.begin(), result_vertices_at_face_18.end(), 24) != result_vertices_at_face_18.end());
      IT_ result_vertices_at_face_19(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 19));
      res.push_back(result_vertices_at_face_19.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_19.begin(), result_vertices_at_face_19.end(), 24) != result_vertices_at_face_19.end());
      res.push_back(std::find(result_vertices_at_face_19.begin(), result_vertices_at_face_19.end(), 12) != result_vertices_at_face_19.end());
      res.push_back(std::find(result_vertices_at_face_19.begin(), result_vertices_at_face_19.end(), 4) != result_vertices_at_face_19.end());
      res.push_back(std::find(result_vertices_at_face_19.begin(), result_vertices_at_face_19.end(), 15) != result_vertices_at_face_19.end());
      IT_ result_vertices_at_face_20(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 20));
      res.push_back(result_vertices_at_face_20.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_20.begin(), result_vertices_at_face_20.end(), 18) != result_vertices_at_face_20.end());
      res.push_back(std::find(result_vertices_at_face_20.begin(), result_vertices_at_face_20.end(), 24) != result_vertices_at_face_20.end());
      res.push_back(std::find(result_vertices_at_face_20.begin(), result_vertices_at_face_20.end(), 15) != result_vertices_at_face_20.end());
      res.push_back(std::find(result_vertices_at_face_20.begin(), result_vertices_at_face_20.end(), 6) != result_vertices_at_face_20.end());
      IT_ result_vertices_at_face_21(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 21));
      res.push_back(result_vertices_at_face_21.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_21.begin(), result_vertices_at_face_21.end(), 9) != result_vertices_at_face_21.end());
      res.push_back(std::find(result_vertices_at_face_21.begin(), result_vertices_at_face_21.end(), 3) != result_vertices_at_face_21.end());
      res.push_back(std::find(result_vertices_at_face_21.begin(), result_vertices_at_face_21.end(), 13) != result_vertices_at_face_21.end());
      res.push_back(std::find(result_vertices_at_face_21.begin(), result_vertices_at_face_21.end(), 25) != result_vertices_at_face_21.end());
      IT_ result_vertices_at_face_22(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 22));
      res.push_back(result_vertices_at_face_22.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_22.begin(), result_vertices_at_face_22.end(), 25) != result_vertices_at_face_22.end());
      res.push_back(std::find(result_vertices_at_face_22.begin(), result_vertices_at_face_22.end(), 13) != result_vertices_at_face_22.end());
      res.push_back(std::find(result_vertices_at_face_22.begin(), result_vertices_at_face_22.end(), 5) != result_vertices_at_face_22.end());
      res.push_back(std::find(result_vertices_at_face_22.begin(), result_vertices_at_face_22.end(), 16) != result_vertices_at_face_22.end());
      IT_ result_vertices_at_face_23(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 23));
      res.push_back(result_vertices_at_face_23.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_23.begin(), result_vertices_at_face_23.end(), 19) != result_vertices_at_face_23.end());
      res.push_back(std::find(result_vertices_at_face_23.begin(), result_vertices_at_face_23.end(), 25) != result_vertices_at_face_23.end());
      res.push_back(std::find(result_vertices_at_face_23.begin(), result_vertices_at_face_23.end(), 16) != result_vertices_at_face_23.end());
      res.push_back(std::find(result_vertices_at_face_23.begin(), result_vertices_at_face_23.end(), 7) != result_vertices_at_face_23.end());
      IT_ result_vertices_at_face_24(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 24));
      res.push_back(result_vertices_at_face_24.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_24.begin(), result_vertices_at_face_24.end(), 20) != result_vertices_at_face_24.end());
      res.push_back(std::find(result_vertices_at_face_24.begin(), result_vertices_at_face_24.end(), 11) != result_vertices_at_face_24.end());
      res.push_back(std::find(result_vertices_at_face_24.begin(), result_vertices_at_face_24.end(), 22) != result_vertices_at_face_24.end());
      res.push_back(std::find(result_vertices_at_face_24.begin(), result_vertices_at_face_24.end(), 26) != result_vertices_at_face_24.end());
      IT_ result_vertices_at_face_25(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 25));
      res.push_back(result_vertices_at_face_25.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_25.begin(), result_vertices_at_face_25.end(), 10) != result_vertices_at_face_25.end());
      res.push_back(std::find(result_vertices_at_face_25.begin(), result_vertices_at_face_25.end(), 20) != result_vertices_at_face_25.end());
      res.push_back(std::find(result_vertices_at_face_25.begin(), result_vertices_at_face_25.end(), 26) != result_vertices_at_face_25.end());
      res.push_back(std::find(result_vertices_at_face_25.begin(), result_vertices_at_face_25.end(), 23) != result_vertices_at_face_25.end());
      IT_ result_vertices_at_face_26(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 26));
      res.push_back(result_vertices_at_face_26.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_26.begin(), result_vertices_at_face_26.end(), 8) != result_vertices_at_face_26.end());
      res.push_back(std::find(result_vertices_at_face_26.begin(), result_vertices_at_face_26.end(), 24) != result_vertices_at_face_26.end());
      res.push_back(std::find(result_vertices_at_face_26.begin(), result_vertices_at_face_26.end(), 26) != result_vertices_at_face_26.end());
      res.push_back(std::find(result_vertices_at_face_26.begin(), result_vertices_at_face_26.end(), 20) != result_vertices_at_face_26.end());
      IT_ result_vertices_at_face_27(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 27));
      res.push_back(result_vertices_at_face_27.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_27.begin(), result_vertices_at_face_27.end(), 20) != result_vertices_at_face_27.end());
      res.push_back(std::find(result_vertices_at_face_27.begin(), result_vertices_at_face_27.end(), 26) != result_vertices_at_face_27.end());
      res.push_back(std::find(result_vertices_at_face_27.begin(), result_vertices_at_face_27.end(), 25) != result_vertices_at_face_27.end());
      res.push_back(std::find(result_vertices_at_face_27.begin(), result_vertices_at_face_27.end(), 9) != result_vertices_at_face_27.end());
      IT_ result_vertices_at_face_28(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 28));
      res.push_back(result_vertices_at_face_28.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_28.begin(), result_vertices_at_face_28.end(), 26) != result_vertices_at_face_28.end());
      res.push_back(std::find(result_vertices_at_face_28.begin(), result_vertices_at_face_28.end(), 22) != result_vertices_at_face_28.end());
      res.push_back(std::find(result_vertices_at_face_28.begin(), result_vertices_at_face_28.end(), 14) != result_vertices_at_face_28.end());
      res.push_back(std::find(result_vertices_at_face_28.begin(), result_vertices_at_face_28.end(), 21) != result_vertices_at_face_28.end());
      IT_ result_vertices_at_face_29(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 29));
      res.push_back(result_vertices_at_face_29.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_29.begin(), result_vertices_at_face_29.end(), 23) != result_vertices_at_face_29.end());
      res.push_back(std::find(result_vertices_at_face_29.begin(), result_vertices_at_face_29.end(), 26) != result_vertices_at_face_29.end());
      res.push_back(std::find(result_vertices_at_face_29.begin(), result_vertices_at_face_29.end(), 21) != result_vertices_at_face_29.end());
      res.push_back(std::find(result_vertices_at_face_29.begin(), result_vertices_at_face_29.end(), 17) != result_vertices_at_face_29.end());
      IT_ result_vertices_at_face_30(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 30));
      res.push_back(result_vertices_at_face_30.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_30.begin(), result_vertices_at_face_30.end(), 24) != result_vertices_at_face_30.end());
      res.push_back(std::find(result_vertices_at_face_30.begin(), result_vertices_at_face_30.end(), 15) != result_vertices_at_face_30.end());
      res.push_back(std::find(result_vertices_at_face_30.begin(), result_vertices_at_face_30.end(), 21) != result_vertices_at_face_30.end());
      res.push_back(std::find(result_vertices_at_face_30.begin(), result_vertices_at_face_30.end(), 26) != result_vertices_at_face_30.end());
      IT_ result_vertices_at_face_31(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 31));
      res.push_back(result_vertices_at_face_31.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_31.begin(), result_vertices_at_face_31.end(), 25) != result_vertices_at_face_31.end());
      res.push_back(std::find(result_vertices_at_face_31.begin(), result_vertices_at_face_31.end(), 26) != result_vertices_at_face_31.end());
      res.push_back(std::find(result_vertices_at_face_31.begin(), result_vertices_at_face_31.end(), 21) != result_vertices_at_face_31.end());
      res.push_back(std::find(result_vertices_at_face_31.begin(), result_vertices_at_face_31.end(), 16) != result_vertices_at_face_31.end());
      IT_ result_vertices_at_face_32(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 32));
      res.push_back(result_vertices_at_face_32.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_32.begin(), result_vertices_at_face_32.end(), 24) != result_vertices_at_face_32.end());
      res.push_back(std::find(result_vertices_at_face_32.begin(), result_vertices_at_face_32.end(), 12) != result_vertices_at_face_32.end());
      res.push_back(std::find(result_vertices_at_face_32.begin(), result_vertices_at_face_32.end(), 22) != result_vertices_at_face_32.end());
      res.push_back(std::find(result_vertices_at_face_32.begin(), result_vertices_at_face_32.end(), 26) != result_vertices_at_face_32.end());
      IT_ result_vertices_at_face_33(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 33));
      res.push_back(result_vertices_at_face_33.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_33.begin(), result_vertices_at_face_33.end(), 26) != result_vertices_at_face_33.end());
      res.push_back(std::find(result_vertices_at_face_33.begin(), result_vertices_at_face_33.end(), 22) != result_vertices_at_face_33.end());
      res.push_back(std::find(result_vertices_at_face_33.begin(), result_vertices_at_face_33.end(), 13) != result_vertices_at_face_33.end());
      res.push_back(std::find(result_vertices_at_face_33.begin(), result_vertices_at_face_33.end(), 25) != result_vertices_at_face_33.end());
      IT_ result_vertices_at_face_34(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 34));
      res.push_back(result_vertices_at_face_34.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_34.begin(), result_vertices_at_face_34.end(), 18) != result_vertices_at_face_34.end());
      res.push_back(std::find(result_vertices_at_face_34.begin(), result_vertices_at_face_34.end(), 24) != result_vertices_at_face_34.end());
      res.push_back(std::find(result_vertices_at_face_34.begin(), result_vertices_at_face_34.end(), 26) != result_vertices_at_face_34.end());
      res.push_back(std::find(result_vertices_at_face_34.begin(), result_vertices_at_face_34.end(), 23) != result_vertices_at_face_34.end());
      IT_ result_vertices_at_face_35(m_fine.get_adjacent_polytopes(pl_face, pl_vertex, 35));
      res.push_back(result_vertices_at_face_35.size() ==  4ul);
      res.push_back(std::find(result_vertices_at_face_35.begin(), result_vertices_at_face_35.end(), 23) != result_vertices_at_face_35.end());
      res.push_back(std::find(result_vertices_at_face_35.begin(), result_vertices_at_face_35.end(), 26) != result_vertices_at_face_35.end());
      res.push_back(std::find(result_vertices_at_face_35.begin(), result_vertices_at_face_35.end(), 25) != result_vertices_at_face_35.end());
      res.push_back(std::find(result_vertices_at_face_35.begin(), result_vertices_at_face_35.end(), 19) != result_vertices_at_face_35.end());

      IT_ result_vertices_at_polyhedron_0(m_fine.get_adjacent_polytopes(pl_polyhedron, pl_vertex, 0));
      res.push_back(result_vertices_at_polyhedron_0.size() ==  8ul);
      res.push_back(std::find(result_vertices_at_polyhedron_0.begin(), result_vertices_at_polyhedron_0.end(), 0) != result_vertices_at_polyhedron_0.end());
      res.push_back(std::find(result_vertices_at_polyhedron_0.begin(), result_vertices_at_polyhedron_0.end(), 8) != result_vertices_at_polyhedron_0.end());
      res.push_back(std::find(result_vertices_at_polyhedron_0.begin(), result_vertices_at_polyhedron_0.end(), 18) != result_vertices_at_polyhedron_0.end());
      res.push_back(std::find(result_vertices_at_polyhedron_0.begin(), result_vertices_at_polyhedron_0.end(), 24) != result_vertices_at_polyhedron_0.end());
      res.push_back(std::find(result_vertices_at_polyhedron_0.begin(), result_vertices_at_polyhedron_0.end(), 10) != result_vertices_at_polyhedron_0.end());
      res.push_back(std::find(result_vertices_at_polyhedron_0.begin(), result_vertices_at_polyhedron_0.end(), 20) != result_vertices_at_polyhedron_0.end());
      res.push_back(std::find(result_vertices_at_polyhedron_0.begin(), result_vertices_at_polyhedron_0.end(), 23) != result_vertices_at_polyhedron_0.end());
      res.push_back(std::find(result_vertices_at_polyhedron_0.begin(), result_vertices_at_polyhedron_0.end(), 26) != result_vertices_at_polyhedron_0.end());
      IT_ result_vertices_at_polyhedron_1(m_fine.get_adjacent_polytopes(pl_polyhedron, pl_vertex, 1));
      res.push_back(result_vertices_at_polyhedron_1.size() ==  8ul);
      res.push_back(std::find(result_vertices_at_polyhedron_1.begin(), result_vertices_at_polyhedron_1.end(), 8) != result_vertices_at_polyhedron_1.end());
      res.push_back(std::find(result_vertices_at_polyhedron_1.begin(), result_vertices_at_polyhedron_1.end(), 1) != result_vertices_at_polyhedron_1.end());
      res.push_back(std::find(result_vertices_at_polyhedron_1.begin(), result_vertices_at_polyhedron_1.end(), 12) != result_vertices_at_polyhedron_1.end());
      res.push_back(std::find(result_vertices_at_polyhedron_1.begin(), result_vertices_at_polyhedron_1.end(), 24) != result_vertices_at_polyhedron_1.end());
      res.push_back(std::find(result_vertices_at_polyhedron_1.begin(), result_vertices_at_polyhedron_1.end(), 20) != result_vertices_at_polyhedron_1.end());
      res.push_back(std::find(result_vertices_at_polyhedron_1.begin(), result_vertices_at_polyhedron_1.end(), 11) != result_vertices_at_polyhedron_1.end());
      res.push_back(std::find(result_vertices_at_polyhedron_1.begin(), result_vertices_at_polyhedron_1.end(), 22) != result_vertices_at_polyhedron_1.end());
      res.push_back(std::find(result_vertices_at_polyhedron_1.begin(), result_vertices_at_polyhedron_1.end(), 26) != result_vertices_at_polyhedron_1.end());
      IT_ result_vertices_at_polyhedron_2(m_fine.get_adjacent_polytopes(pl_polyhedron, pl_vertex, 2));
      res.push_back(result_vertices_at_polyhedron_2.size() ==  8ul);
      res.push_back(std::find(result_vertices_at_polyhedron_2.begin(), result_vertices_at_polyhedron_2.end(), 10) != result_vertices_at_polyhedron_2.end());
      res.push_back(std::find(result_vertices_at_polyhedron_2.begin(), result_vertices_at_polyhedron_2.end(), 20) != result_vertices_at_polyhedron_2.end());
      res.push_back(std::find(result_vertices_at_polyhedron_2.begin(), result_vertices_at_polyhedron_2.end(), 26) != result_vertices_at_polyhedron_2.end());
      res.push_back(std::find(result_vertices_at_polyhedron_2.begin(), result_vertices_at_polyhedron_2.end(), 23) != result_vertices_at_polyhedron_2.end());
      res.push_back(std::find(result_vertices_at_polyhedron_2.begin(), result_vertices_at_polyhedron_2.end(), 2) != result_vertices_at_polyhedron_2.end());
      res.push_back(std::find(result_vertices_at_polyhedron_2.begin(), result_vertices_at_polyhedron_2.end(), 9) != result_vertices_at_polyhedron_2.end());
      res.push_back(std::find(result_vertices_at_polyhedron_2.begin(), result_vertices_at_polyhedron_2.end(), 25) != result_vertices_at_polyhedron_2.end());
      res.push_back(std::find(result_vertices_at_polyhedron_2.begin(), result_vertices_at_polyhedron_2.end(), 19) != result_vertices_at_polyhedron_2.end());
      IT_ result_vertices_at_polyhedron_3(m_fine.get_adjacent_polytopes(pl_polyhedron, pl_vertex, 3));
      res.push_back(result_vertices_at_polyhedron_3.size() ==  8ul);
      res.push_back(std::find(result_vertices_at_polyhedron_3.begin(), result_vertices_at_polyhedron_3.end(), 20) != result_vertices_at_polyhedron_3.end());
      res.push_back(std::find(result_vertices_at_polyhedron_3.begin(), result_vertices_at_polyhedron_3.end(), 11) != result_vertices_at_polyhedron_3.end());
      res.push_back(std::find(result_vertices_at_polyhedron_3.begin(), result_vertices_at_polyhedron_3.end(), 22) != result_vertices_at_polyhedron_3.end());
      res.push_back(std::find(result_vertices_at_polyhedron_3.begin(), result_vertices_at_polyhedron_3.end(), 26) != result_vertices_at_polyhedron_3.end());
      res.push_back(std::find(result_vertices_at_polyhedron_3.begin(), result_vertices_at_polyhedron_3.end(), 9) != result_vertices_at_polyhedron_3.end());
      res.push_back(std::find(result_vertices_at_polyhedron_3.begin(), result_vertices_at_polyhedron_3.end(), 3) != result_vertices_at_polyhedron_3.end());
      res.push_back(std::find(result_vertices_at_polyhedron_3.begin(), result_vertices_at_polyhedron_3.end(), 13) != result_vertices_at_polyhedron_3.end());
      res.push_back(std::find(result_vertices_at_polyhedron_3.begin(), result_vertices_at_polyhedron_3.end(), 25) != result_vertices_at_polyhedron_3.end());
      IT_ result_vertices_at_polyhedron_4(m_fine.get_adjacent_polytopes(pl_polyhedron, pl_vertex, 4));
      res.push_back(result_vertices_at_polyhedron_4.size() ==  8ul);
      res.push_back(std::find(result_vertices_at_polyhedron_4.begin(), result_vertices_at_polyhedron_4.end(), 24) != result_vertices_at_polyhedron_4.end());
      res.push_back(std::find(result_vertices_at_polyhedron_4.begin(), result_vertices_at_polyhedron_4.end(), 12) != result_vertices_at_polyhedron_4.end());
      res.push_back(std::find(result_vertices_at_polyhedron_4.begin(), result_vertices_at_polyhedron_4.end(), 4) != result_vertices_at_polyhedron_4.end());
      res.push_back(std::find(result_vertices_at_polyhedron_4.begin(), result_vertices_at_polyhedron_4.end(), 15) != result_vertices_at_polyhedron_4.end());
      res.push_back(std::find(result_vertices_at_polyhedron_4.begin(), result_vertices_at_polyhedron_4.end(), 26) != result_vertices_at_polyhedron_4.end());
      res.push_back(std::find(result_vertices_at_polyhedron_4.begin(), result_vertices_at_polyhedron_4.end(), 22) != result_vertices_at_polyhedron_4.end());
      res.push_back(std::find(result_vertices_at_polyhedron_4.begin(), result_vertices_at_polyhedron_4.end(), 14) != result_vertices_at_polyhedron_4.end());
      res.push_back(std::find(result_vertices_at_polyhedron_4.begin(), result_vertices_at_polyhedron_4.end(), 21) != result_vertices_at_polyhedron_4.end());
      IT_ result_vertices_at_polyhedron_5(m_fine.get_adjacent_polytopes(pl_polyhedron, pl_vertex, 5));
      res.push_back(result_vertices_at_polyhedron_5.size() ==  8ul);
      res.push_back(std::find(result_vertices_at_polyhedron_5.begin(), result_vertices_at_polyhedron_5.end(), 26) != result_vertices_at_polyhedron_5.end());
      res.push_back(std::find(result_vertices_at_polyhedron_5.begin(), result_vertices_at_polyhedron_5.end(), 22) != result_vertices_at_polyhedron_5.end());
      res.push_back(std::find(result_vertices_at_polyhedron_5.begin(), result_vertices_at_polyhedron_5.end(), 14) != result_vertices_at_polyhedron_5.end());
      res.push_back(std::find(result_vertices_at_polyhedron_5.begin(), result_vertices_at_polyhedron_5.end(), 21) != result_vertices_at_polyhedron_5.end());
      res.push_back(std::find(result_vertices_at_polyhedron_5.begin(), result_vertices_at_polyhedron_5.end(), 25) != result_vertices_at_polyhedron_5.end());
      res.push_back(std::find(result_vertices_at_polyhedron_5.begin(), result_vertices_at_polyhedron_5.end(), 13) != result_vertices_at_polyhedron_5.end());
      res.push_back(std::find(result_vertices_at_polyhedron_5.begin(), result_vertices_at_polyhedron_5.end(), 5) != result_vertices_at_polyhedron_5.end());
      res.push_back(std::find(result_vertices_at_polyhedron_5.begin(), result_vertices_at_polyhedron_5.end(), 16) != result_vertices_at_polyhedron_5.end());
      IT_ result_vertices_at_polyhedron_6(m_fine.get_adjacent_polytopes(pl_polyhedron, pl_vertex, 6));
      res.push_back(result_vertices_at_polyhedron_6.size() ==  8ul);
      res.push_back(std::find(result_vertices_at_polyhedron_6.begin(), result_vertices_at_polyhedron_6.end(), 18) != result_vertices_at_polyhedron_6.end());
      res.push_back(std::find(result_vertices_at_polyhedron_6.begin(), result_vertices_at_polyhedron_6.end(), 24) != result_vertices_at_polyhedron_6.end());
      res.push_back(std::find(result_vertices_at_polyhedron_6.begin(), result_vertices_at_polyhedron_6.end(), 15) != result_vertices_at_polyhedron_6.end());
      res.push_back(std::find(result_vertices_at_polyhedron_6.begin(), result_vertices_at_polyhedron_6.end(), 6) != result_vertices_at_polyhedron_6.end());
      res.push_back(std::find(result_vertices_at_polyhedron_6.begin(), result_vertices_at_polyhedron_6.end(), 23) != result_vertices_at_polyhedron_6.end());
      res.push_back(std::find(result_vertices_at_polyhedron_6.begin(), result_vertices_at_polyhedron_6.end(), 26) != result_vertices_at_polyhedron_6.end());
      res.push_back(std::find(result_vertices_at_polyhedron_6.begin(), result_vertices_at_polyhedron_6.end(), 21) != result_vertices_at_polyhedron_6.end());
      res.push_back(std::find(result_vertices_at_polyhedron_6.begin(), result_vertices_at_polyhedron_6.end(), 17) != result_vertices_at_polyhedron_6.end());
      IT_ result_vertices_at_polyhedron_7(m_fine.get_adjacent_polytopes(pl_polyhedron, pl_vertex, 7));
      res.push_back(result_vertices_at_polyhedron_7.size() ==  8ul);
      res.push_back(std::find(result_vertices_at_polyhedron_7.begin(), result_vertices_at_polyhedron_7.end(), 23) != result_vertices_at_polyhedron_7.end());
      res.push_back(std::find(result_vertices_at_polyhedron_7.begin(), result_vertices_at_polyhedron_7.end(), 26) != result_vertices_at_polyhedron_7.end());
      res.push_back(std::find(result_vertices_at_polyhedron_7.begin(), result_vertices_at_polyhedron_7.end(), 21) != result_vertices_at_polyhedron_7.end());
      res.push_back(std::find(result_vertices_at_polyhedron_7.begin(), result_vertices_at_polyhedron_7.end(), 17) != result_vertices_at_polyhedron_7.end());
      res.push_back(std::find(result_vertices_at_polyhedron_7.begin(), result_vertices_at_polyhedron_7.end(), 19) != result_vertices_at_polyhedron_7.end());
      res.push_back(std::find(result_vertices_at_polyhedron_7.begin(), result_vertices_at_polyhedron_7.end(), 25) != result_vertices_at_polyhedron_7.end());
      res.push_back(std::find(result_vertices_at_polyhedron_7.begin(), result_vertices_at_polyhedron_7.end(), 16) != result_vertices_at_polyhedron_7.end());
      res.push_back(std::find(result_vertices_at_polyhedron_7.begin(), result_vertices_at_polyhedron_7.end(), 7) != result_vertices_at_polyhedron_7.end());


      IT_ result_polyhedron_at_vertex_0(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 0));
      res.push_back(result_polyhedron_at_vertex_0.size() ==  1ul);
      res.push_back(std::find(result_polyhedron_at_vertex_0.begin(), result_polyhedron_at_vertex_0.end(), 0) != result_polyhedron_at_vertex_0.end());

      IT_ result_polyhedron_at_vertex_1(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 1));
      res.push_back(result_polyhedron_at_vertex_1.size() ==  1ul);
      res.push_back(std::find(result_polyhedron_at_vertex_1.begin(), result_polyhedron_at_vertex_1.end(), 1) != result_polyhedron_at_vertex_1.end());

      IT_ result_polyhedron_at_vertex_2(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 2));
      res.push_back(result_polyhedron_at_vertex_2.size() ==  1ul);
      res.push_back(std::find(result_polyhedron_at_vertex_2.begin(), result_polyhedron_at_vertex_2.end(), 2) != result_polyhedron_at_vertex_2.end());

      IT_ result_polyhedron_at_vertex_3(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 3));
      res.push_back(result_polyhedron_at_vertex_3.size() ==  1ul);
      res.push_back(std::find(result_polyhedron_at_vertex_3.begin(), result_polyhedron_at_vertex_3.end(), 3) != result_polyhedron_at_vertex_3.end());

      IT_ result_polyhedron_at_vertex_4(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 4));
      res.push_back(result_polyhedron_at_vertex_4.size() ==  1ul);
      res.push_back(std::find(result_polyhedron_at_vertex_4.begin(), result_polyhedron_at_vertex_4.end(), 4) != result_polyhedron_at_vertex_4.end());

      IT_ result_polyhedron_at_vertex_5(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 5));
      res.push_back(result_polyhedron_at_vertex_5.size() ==  1ul);
      res.push_back(std::find(result_polyhedron_at_vertex_5.begin(), result_polyhedron_at_vertex_5.end(), 5) != result_polyhedron_at_vertex_5.end());

      IT_ result_polyhedron_at_vertex_6(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 6));
      res.push_back(result_polyhedron_at_vertex_6.size() ==  1ul);
      res.push_back(std::find(result_polyhedron_at_vertex_6.begin(), result_polyhedron_at_vertex_6.end(), 6) != result_polyhedron_at_vertex_6.end());

      IT_ result_polyhedron_at_vertex_7(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 7));
      res.push_back(result_polyhedron_at_vertex_7.size() ==  1ul);
      res.push_back(std::find(result_polyhedron_at_vertex_7.begin(), result_polyhedron_at_vertex_7.end(), 7) != result_polyhedron_at_vertex_7.end());

      IT_ result_polyhedron_at_vertex_8(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 8));
      res.push_back(result_polyhedron_at_vertex_8.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_8.begin(), result_polyhedron_at_vertex_8.end(), 0) != result_polyhedron_at_vertex_8.end());
      res.push_back(std::find(result_polyhedron_at_vertex_8.begin(), result_polyhedron_at_vertex_8.end(), 1) != result_polyhedron_at_vertex_8.end());

      IT_ result_polyhedron_at_vertex_9(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 9));
      res.push_back(result_polyhedron_at_vertex_9.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_9.begin(), result_polyhedron_at_vertex_9.end(), 2) != result_polyhedron_at_vertex_9.end());
      res.push_back(std::find(result_polyhedron_at_vertex_9.begin(), result_polyhedron_at_vertex_9.end(), 3) != result_polyhedron_at_vertex_9.end());

      IT_ result_polyhedron_at_vertex_10(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 10));
      res.push_back(result_polyhedron_at_vertex_10.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_10.begin(), result_polyhedron_at_vertex_10.end(), 0) != result_polyhedron_at_vertex_10.end());
      res.push_back(std::find(result_polyhedron_at_vertex_10.begin(), result_polyhedron_at_vertex_10.end(), 2) != result_polyhedron_at_vertex_10.end());

      IT_ result_polyhedron_at_vertex_11(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 11));
      res.push_back(result_polyhedron_at_vertex_11.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_11.begin(), result_polyhedron_at_vertex_11.end(), 1) != result_polyhedron_at_vertex_11.end());
      res.push_back(std::find(result_polyhedron_at_vertex_11.begin(), result_polyhedron_at_vertex_11.end(), 3) != result_polyhedron_at_vertex_11.end());

      IT_ result_polyhedron_at_vertex_12(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 12));
      res.push_back(result_polyhedron_at_vertex_12.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_12.begin(), result_polyhedron_at_vertex_12.end(), 1) != result_polyhedron_at_vertex_12.end());
      res.push_back(std::find(result_polyhedron_at_vertex_12.begin(), result_polyhedron_at_vertex_12.end(), 4) != result_polyhedron_at_vertex_12.end());

      IT_ result_polyhedron_at_vertex_13(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 13));
      res.push_back(result_polyhedron_at_vertex_13.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_13.begin(), result_polyhedron_at_vertex_13.end(), 3) != result_polyhedron_at_vertex_13.end());
      res.push_back(std::find(result_polyhedron_at_vertex_13.begin(), result_polyhedron_at_vertex_13.end(), 5) != result_polyhedron_at_vertex_13.end());

      IT_ result_polyhedron_at_vertex_14(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 14));
      res.push_back(result_polyhedron_at_vertex_14.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_14.begin(), result_polyhedron_at_vertex_14.end(), 4) != result_polyhedron_at_vertex_14.end());
      res.push_back(std::find(result_polyhedron_at_vertex_14.begin(), result_polyhedron_at_vertex_14.end(), 5) != result_polyhedron_at_vertex_14.end());

      IT_ result_polyhedron_at_vertex_15(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 15));
      res.push_back(result_polyhedron_at_vertex_15.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_15.begin(), result_polyhedron_at_vertex_15.end(), 4) != result_polyhedron_at_vertex_15.end());
      res.push_back(std::find(result_polyhedron_at_vertex_15.begin(), result_polyhedron_at_vertex_15.end(), 6) != result_polyhedron_at_vertex_15.end());

      IT_ result_polyhedron_at_vertex_16(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 16));
      res.push_back(result_polyhedron_at_vertex_16.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_16.begin(), result_polyhedron_at_vertex_16.end(), 5) != result_polyhedron_at_vertex_16.end());
      res.push_back(std::find(result_polyhedron_at_vertex_16.begin(), result_polyhedron_at_vertex_16.end(), 7) != result_polyhedron_at_vertex_16.end());

      IT_ result_polyhedron_at_vertex_17(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 17));
      res.push_back(result_polyhedron_at_vertex_17.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_17.begin(), result_polyhedron_at_vertex_17.end(), 6) != result_polyhedron_at_vertex_17.end());
      res.push_back(std::find(result_polyhedron_at_vertex_17.begin(), result_polyhedron_at_vertex_17.end(), 7) != result_polyhedron_at_vertex_17.end());

      IT_ result_polyhedron_at_vertex_18(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 18));
      res.push_back(result_polyhedron_at_vertex_18.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_18.begin(), result_polyhedron_at_vertex_18.end(), 0) != result_polyhedron_at_vertex_18.end());
      res.push_back(std::find(result_polyhedron_at_vertex_18.begin(), result_polyhedron_at_vertex_18.end(), 6) != result_polyhedron_at_vertex_18.end());

      IT_ result_polyhedron_at_vertex_19(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 19));
      res.push_back(result_polyhedron_at_vertex_19.size() ==  2ul);
      res.push_back(std::find(result_polyhedron_at_vertex_19.begin(), result_polyhedron_at_vertex_19.end(), 2) != result_polyhedron_at_vertex_19.end());
      res.push_back(std::find(result_polyhedron_at_vertex_19.begin(), result_polyhedron_at_vertex_19.end(), 7) != result_polyhedron_at_vertex_19.end());

      IT_ result_polyhedron_at_vertex_20(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 20));
      res.push_back(result_polyhedron_at_vertex_20.size() ==  4ul);
      res.push_back(std::find(result_polyhedron_at_vertex_20.begin(), result_polyhedron_at_vertex_20.end(), 0) != result_polyhedron_at_vertex_20.end());
      res.push_back(std::find(result_polyhedron_at_vertex_20.begin(), result_polyhedron_at_vertex_20.end(), 1) != result_polyhedron_at_vertex_20.end());
      res.push_back(std::find(result_polyhedron_at_vertex_20.begin(), result_polyhedron_at_vertex_20.end(), 2) != result_polyhedron_at_vertex_20.end());
      res.push_back(std::find(result_polyhedron_at_vertex_20.begin(), result_polyhedron_at_vertex_20.end(), 3) != result_polyhedron_at_vertex_20.end());

      IT_ result_polyhedron_at_vertex_21(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 21));
      res.push_back(result_polyhedron_at_vertex_21.size() ==  4ul);
      res.push_back(std::find(result_polyhedron_at_vertex_21.begin(), result_polyhedron_at_vertex_21.end(), 4) != result_polyhedron_at_vertex_21.end());
      res.push_back(std::find(result_polyhedron_at_vertex_21.begin(), result_polyhedron_at_vertex_21.end(), 5) != result_polyhedron_at_vertex_21.end());
      res.push_back(std::find(result_polyhedron_at_vertex_21.begin(), result_polyhedron_at_vertex_21.end(), 6) != result_polyhedron_at_vertex_21.end());
      res.push_back(std::find(result_polyhedron_at_vertex_21.begin(), result_polyhedron_at_vertex_21.end(), 7) != result_polyhedron_at_vertex_21.end());

      IT_ result_polyhedron_at_vertex_22(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 22));
      res.push_back(result_polyhedron_at_vertex_22.size() ==  4ul);
      res.push_back(std::find(result_polyhedron_at_vertex_22.begin(), result_polyhedron_at_vertex_22.end(), 1) != result_polyhedron_at_vertex_22.end());
      res.push_back(std::find(result_polyhedron_at_vertex_22.begin(), result_polyhedron_at_vertex_22.end(), 3) != result_polyhedron_at_vertex_22.end());
      res.push_back(std::find(result_polyhedron_at_vertex_22.begin(), result_polyhedron_at_vertex_22.end(), 4) != result_polyhedron_at_vertex_22.end());
      res.push_back(std::find(result_polyhedron_at_vertex_22.begin(), result_polyhedron_at_vertex_22.end(), 5) != result_polyhedron_at_vertex_22.end());

      IT_ result_polyhedron_at_vertex_23(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 23));
      res.push_back(result_polyhedron_at_vertex_23.size() ==  4ul);
      res.push_back(std::find(result_polyhedron_at_vertex_23.begin(), result_polyhedron_at_vertex_23.end(), 0) != result_polyhedron_at_vertex_23.end());
      res.push_back(std::find(result_polyhedron_at_vertex_23.begin(), result_polyhedron_at_vertex_23.end(), 2) != result_polyhedron_at_vertex_23.end());
      res.push_back(std::find(result_polyhedron_at_vertex_23.begin(), result_polyhedron_at_vertex_23.end(), 6) != result_polyhedron_at_vertex_23.end());
      res.push_back(std::find(result_polyhedron_at_vertex_23.begin(), result_polyhedron_at_vertex_23.end(), 7) != result_polyhedron_at_vertex_23.end());

      IT_ result_polyhedron_at_vertex_24(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 24));
      res.push_back(result_polyhedron_at_vertex_24.size() ==  4ul);
      res.push_back(std::find(result_polyhedron_at_vertex_24.begin(), result_polyhedron_at_vertex_24.end(), 0) != result_polyhedron_at_vertex_24.end());
      res.push_back(std::find(result_polyhedron_at_vertex_24.begin(), result_polyhedron_at_vertex_24.end(), 1) != result_polyhedron_at_vertex_24.end());
      res.push_back(std::find(result_polyhedron_at_vertex_24.begin(), result_polyhedron_at_vertex_24.end(), 4) != result_polyhedron_at_vertex_24.end());
      res.push_back(std::find(result_polyhedron_at_vertex_24.begin(), result_polyhedron_at_vertex_24.end(), 6) != result_polyhedron_at_vertex_24.end());

      IT_ result_polyhedron_at_vertex_25(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 25));
      res.push_back(result_polyhedron_at_vertex_25.size() ==  4ul);
      res.push_back(std::find(result_polyhedron_at_vertex_25.begin(), result_polyhedron_at_vertex_25.end(), 2) != result_polyhedron_at_vertex_25.end());
      res.push_back(std::find(result_polyhedron_at_vertex_25.begin(), result_polyhedron_at_vertex_25.end(), 3) != result_polyhedron_at_vertex_25.end());
      res.push_back(std::find(result_polyhedron_at_vertex_25.begin(), result_polyhedron_at_vertex_25.end(), 5) != result_polyhedron_at_vertex_25.end());
      res.push_back(std::find(result_polyhedron_at_vertex_25.begin(), result_polyhedron_at_vertex_25.end(), 7) != result_polyhedron_at_vertex_25.end());

      IT_ result_polyhedron_at_vertex_26(m_fine.get_adjacent_polytopes(pl_vertex, pl_polyhedron, 26));
      res.push_back(result_polyhedron_at_vertex_26.size() ==  8ul);
      res.push_back(std::find(result_polyhedron_at_vertex_26.begin(), result_polyhedron_at_vertex_26.end(), 0) != result_polyhedron_at_vertex_26.end());
      res.push_back(std::find(result_polyhedron_at_vertex_26.begin(), result_polyhedron_at_vertex_26.end(), 1) != result_polyhedron_at_vertex_26.end());
      res.push_back(std::find(result_polyhedron_at_vertex_26.begin(), result_polyhedron_at_vertex_26.end(), 2) != result_polyhedron_at_vertex_26.end());
      res.push_back(std::find(result_polyhedron_at_vertex_26.begin(), result_polyhedron_at_vertex_26.end(), 3) != result_polyhedron_at_vertex_26.end());
      res.push_back(std::find(result_polyhedron_at_vertex_26.begin(), result_polyhedron_at_vertex_26.end(), 4) != result_polyhedron_at_vertex_26.end());
      res.push_back(std::find(result_polyhedron_at_vertex_26.begin(), result_polyhedron_at_vertex_26.end(), 5) != result_polyhedron_at_vertex_26.end());
      res.push_back(std::find(result_polyhedron_at_vertex_26.begin(), result_polyhedron_at_vertex_26.end(), 6) != result_polyhedron_at_vertex_26.end());
      res.push_back(std::find(result_polyhedron_at_vertex_26.begin(), result_polyhedron_at_vertex_26.end(), 7) != result_polyhedron_at_vertex_26.end());

      res.push_back(halos.at(0)->get_elements().size() ==  IndexType_(8));
      res.push_back(halos.at(0)->get_elements().at(0) ==  IndexType_(0));
      res.push_back(halos.at(0)->get_elements().at(1) ==  IndexType_(7));
      res.push_back(halos.at(0)->get_elements().at(2) ==  IndexType_(6));
      res.push_back(halos.at(0)->get_elements().at(3) ==  IndexType_(5));
      res.push_back(halos.at(0)->get_elements().at(4) ==  IndexType_(4));
      res.push_back(halos.at(0)->get_elements().at(5) ==  IndexType_(3));
      res.push_back(halos.at(0)->get_elements().at(6) ==  IndexType_(2));
      res.push_back(halos.at(0)->get_elements().at(7) ==  IndexType_(1));

      res.push_back(halos.at(1)->get_elements().size() ==  IndexType_(4));
      res.push_back(halos.at(1)->get_elements().at(0) ==  IndexType_(3));
      res.push_back(halos.at(1)->get_elements().at(1) ==  IndexType_(17));
      res.push_back(halos.at(1)->get_elements().at(2) ==  IndexType_(16));
      res.push_back(halos.at(1)->get_elements().at(3) ==  IndexType_(15));

      res.push_back(halos.at(2)->get_elements().size() ==  IndexType_(2));
      res.push_back(halos.at(2)->get_elements().at(0) ==  IndexType_(3));
      res.push_back(halos.at(2)->get_elements().at(1) ==  IndexType_(15));

      res.push_back(attrs.at(0).get_data().size() ==  27ul);
      res.push_back(attrs.at(1).get_data().size() ==  27ul);
      res.push_back(attrs.at(2).get_data().size() ==  27ul);

      res.push_back(std::abs(attrs.at(0).get_data().at(0) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(0) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(0) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(1) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(1) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(1) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(2) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(2) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(2) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(3) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(3) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(3) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(4) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(4) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(4) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(5) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(5) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(5) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(6) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(6) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(6) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(7) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(7) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(7) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(8) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(8) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(8) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(9) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(9) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(9) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(10) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(10) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(10) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(11) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(11) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(11) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(12) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(12) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(12) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(13) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(13) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(13) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(14) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(14) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(14) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(15) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(15) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(15) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(16) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(16) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(16) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(17) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(17) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(17) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(18) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(18) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(18) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(19) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(19) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(19) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(20) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(20) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(20) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(21) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(21) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(21) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(22) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(22) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(22) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(23) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(23) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(23) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(24) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(24) -  double(0.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(24) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(25) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(25) -  double(1.0)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(25) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(0).get_data().at(26) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(1).get_data().at(26) -  double(0.5)) <  std::numeric_limits<double>::epsilon());
      res.push_back(std::abs(attrs.at(2).get_data().at(26) -  double(0.5)) <  std::numeric_limits<double>::epsilon());

      //testing all find calls simultaneously in order to not overencumber the function inliner
      for(auto res_i : res)
      {
        TEST_CHECK(res_i);
      }
    }
};
RefinementTest3D<Mem::Main, Index, Algo::Generic, std::vector, std::vector<Index> > ref_test3_cpu_v_v("std::vector, std::vector");
RefinementTest3D<Mem::Main, Index, Algo::Generic, std::vector, std::deque<Index> > ref_test3_cpu_v_d("std::vector, std::deque");
RefinementTest3D<Mem::Main, Index, Algo::Generic, std::deque, std::vector<Index> > ref_test3_cpu_d_v("std::deque, std::vector");
RefinementTest3D<Mem::Main, Index, Algo::Generic, std::deque, std::deque<Index> > ref_test3_cpu_d_d("std::deque, std::deque");
