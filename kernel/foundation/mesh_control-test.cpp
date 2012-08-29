#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/mesh_control.hpp>
#include<kernel/foundation/dense_data_wrapper.hpp>
#include<kernel/lafem/dense_vector.hpp>
#include<kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;
using namespace FEAST::Geometry;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MeshControlTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    MeshControlTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("MeshControlTest<" + tag + ">")
    {
    }

    void run() const
    {
      /*  2    3
       *  *-1--*
       *  2    |
       *  |    3
       *  *--0-*
       *  0    1
       */

      //creating foundation mesh
      Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > m(0);
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

      typedef ConformalMesh<ConformalMeshPolicy<Shape::Hypercube<dim_2D> > > confmeshtype_;

      IndexType_* size_set(new IndexType_[3]);
      MeshControl<dim_2D>::fill(m, size_set);

      TEST_CHECK_EQUAL(size_set[0], 4);
      TEST_CHECK_EQUAL(size_set[1], 4);
      TEST_CHECK_EQUAL(size_set[2], 1);

      confmeshtype_ confmesh(size_set);
      MeshControl<dim_2D>::fill(m, confmesh);

      typename confmeshtype_::template IndexSet<1,0>::Type& geo_vertex_at_edge(confmesh.template get_index_set<1,0>());
      typename confmeshtype_::template IndexSet<2,0>::Type& geo_vertex_at_face(confmesh.template get_index_set<2,0>());

      TEST_CHECK_EQUAL(geo_vertex_at_edge[0][0], 0);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[0][1], 1);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[1][0], 2);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[1][1], 3);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[2][0], 0);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[2][1], 2);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[3][0], 1);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[3][1], 3);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][0], 0);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][1], 1);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][2], 2);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][3], 3);

      delete[] size_set;
    }
};
MeshControlTest<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > meshcontrol_testvv("std::vector, std::vector");
MeshControlTest<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > meshcontrol_testvd("std::vector, std::deque");
MeshControlTest<Archs::None, unsigned long, std::deque, std::deque<unsigned long> > meshcontrol_testdd("std::deque, std::deque");
MeshControlTest<Archs::None, unsigned long, std::deque, std::vector<unsigned long> > meshcontrol_testdv("std::deque, std::vector");
