#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/attribute.hpp>
#include<kernel/foundation/sub_mesh.hpp>
#include<kernel/foundation/halo.hpp>
#include<kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

using namespace FEAST::Foundation;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class SubMeshTest :
  public TaggedTest<Tag_, IndexType_>
{
  public:
    SubMeshTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("SubMeshTest<" + tag + ">")
    {
    }

    void run() const override
    {
      /*(0,1) (1,1)
       *  *----*
       *  |    |
       *  |    |
       *  *----*
       *(0,0) (1,0)
       */

      //create attributes for vertex coords
      std::vector<std::shared_ptr<Foundation::AttributeBase<OT_> > > attrs;
      attrs.push_back(std::shared_ptr<Foundation::AttributeBase<OT_> >(new Foundation::Attribute<double, OT_>)); //vertex x-coords
      attrs.push_back(std::shared_ptr<Foundation::AttributeBase<OT_> >(new Foundation::Attribute<double, OT_>)); //vertex y-coords

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(0));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(1));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(0));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(1));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(1));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(1));

      /*  2    3
       *  *-1--*
       *  2    |
       *  |    3
       *  *--0-*
       *  0    1
       */

      //creating foundation mesh
      Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_> > m(0);

      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);

      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);

      m.add_polytope(Foundation::pl_face);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 0, 0);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 0, 2);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 0, 0);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 1, 0);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 1, 3);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 1, 0);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 2, 1);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 2, 2);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 2, 0);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 3, 1);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 3, 3);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 3, 0);

      std::shared_ptr<HaloBase<Mesh<Dim2D, Topology<IndexType_, OT_, IT_> > > > h(new Halo<0, PLEdge, Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_> > >(m));
      h->push_back(3);

      ///SUBDIM
      SubMesh<Dim2D,
              Topology<IndexType_, OT_, IT_>,
              std::vector> sm(h);

      TEST_CHECK_EQUAL(sm.get_topologies().at(ipi_vertex_edge).size(), 2ul);
      TEST_CHECK_EQUAL(sm.get_topologies().at(ipi_edge_vertex).size(), 1ul);
      TEST_CHECK_EQUAL(sm.get_map().size(), 1ul);
      TEST_CHECK_EQUAL(sm.get_map().at(0), 3ul);

      TEST_CHECK_EQUAL(sm.get_topologies().at(ipi_vertex_edge).at(0).at(0), 0ul);
      TEST_CHECK_EQUAL(sm.get_topologies().at(ipi_vertex_edge).at(1).at(0), 0ul);
      TEST_CHECK_EQUAL(sm.get_topologies().at(ipi_edge_vertex).at(0).at(0), 0ul);
      TEST_CHECK_EQUAL(sm.get_topologies().at(ipi_edge_vertex).at(0).at(1), 1ul);

      ///DIM
      std::shared_ptr<HaloBase<Mesh<Dim2D, Topology<IndexType_, OT_, IT_> > > > h_face(new Foundation::Halo<0, PLFace, Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_> > >(m));
      h_face->push_back(0);

      SubMesh<Dim2D,
              Topology<IndexType_, OT_, IT_>,
              std::vector> sm_face(h_face);

      TEST_CHECK_EQUAL(sm_face.get_topologies().at(ipi_vertex_edge).size(), 4ul);
      TEST_CHECK_EQUAL(sm_face.get_topologies().at(ipi_edge_vertex).size(), 4ul);
      TEST_CHECK_EQUAL(sm_face.get_topologies().at(ipi_face_vertex).size(), 1ul);

      TEST_CHECK_EQUAL(sm_face.get_topologies().at(ipi_face_vertex).at(0).size(), 4ul);

      auto adj_f_v(sm_face.get_adjacent_polytopes(pl_face, pl_vertex, 0));
      TEST_CHECK_EQUAL(adj_f_v.size(), 4ul);
      TEST_CHECK_EQUAL(adj_f_v.at(0), 0ul);
      TEST_CHECK_EQUAL(adj_f_v.at(1), 1ul);
      TEST_CHECK_EQUAL(adj_f_v.at(2), 2ul);
      TEST_CHECK_EQUAL(adj_f_v.at(3), 3ul);

      TEST_CHECK_EQUAL(sm_face.get_map().size(), 1ul);
      TEST_CHECK_EQUAL(sm_face.get_map().at(0), 0ul);
    }
};
SubMeshTest<Mem::Main, Index, std::vector, std::vector<Index> > submesh_test_fginter_cpu_v_v("std::vector, std::vector");
SubMeshTest<Mem::Main, Index, std::vector, std::deque<Index> > submesh_test_fginter_cpu_v_d("std::vector, std::deque");
SubMeshTest<Mem::Main, Index, std::deque, std::vector<Index> > submesh_test_fginter_cpu_d_v("std::deque, std::vector");
SubMeshTest<Mem::Main, Index, std::deque, std::deque<Index> > submesh_test_fginter_cpu_d_d("std::deque, std::deque");
