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
class MeshControlTest1D:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    MeshControlTest1D(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("MeshControlTest1D<" + tag + ">")
    {
    }

    void run() const
    {
      /* (0)  (1)
       *  *----*
       */

      //create attributes for vertex coords
      std::vector<std::shared_ptr<Foundation::AttributeBase<OT_> > > attrs;
      attrs.push_back(std::shared_ptr<Foundation::AttributeBase<OT_> >(new Foundation::Attribute<double, OT_>)); //vertex x-coords

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(1)));

      /*
       *  *--0-*
       *  0    1
       */

      //creating foundation mesh
      Foundation::Mesh<Foundation::rnt_1D, Foundation::Topology<IndexType_, OT_, IT_> > m(0, &attrs);
      Foundation::MeshAttributeRegistration::execute(m, Foundation::pl_vertex);

      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);

      m.add_polytope(Foundation::pl_edge);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 0, 0);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 1, 0);

      typedef ConformalMesh<ConformalMeshPolicy<Shape::Hypercube<1> > > confmeshtype_;

      IndexType_* size_set(new IndexType_[2]);
      MeshControl<dim_1D>::fill_sizes(m, size_set);

      TEST_CHECK_EQUAL(size_set[0], 2);
      TEST_CHECK_EQUAL(size_set[1], 1);

      confmeshtype_ confmesh(size_set);
      MeshControl<dim_1D>::fill_adjacencies(m, confmesh);

      typename confmeshtype_::template IndexSet<1,0>::Type& geo_vertex_at_edge(confmesh.template get_index_set<1,0>());
      TEST_CHECK_EQUAL(geo_vertex_at_edge[0][0], 0u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[0][1], 1u);

      MeshControl<dim_1D>::fill_vertex_sets(m, confmesh, *((Attribute<double, OT_>*)(attrs.at(0).get())));
      typename confmeshtype_::VertexSetType& vertex_coord_tuples(confmesh.get_vertex_set());
      TEST_CHECK_EQUAL(vertex_coord_tuples[0][0], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[1][0], double(1));

      delete[] size_set;
    }
};
MeshControlTest1D<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > meshcontrol1d_testvv("std::vector, std::vector");
MeshControlTest1D<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > meshcontrol1d_testvd("std::vector, std::deque");

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MeshControlTest2D:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    MeshControlTest2D(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("MeshControlTest2D<" + tag + ">")
    {
    }

    void run() const
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

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(0)));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(1)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(0)));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(1)));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(1)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(1)));

      /*  2    3
       *  *-1--*
       *  2    |
       *  |    3
       *  *--0-*
       *  0    1
       */

      //creating foundation mesh
      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m(0, &attrs);
      Foundation::MeshAttributeRegistration::execute(m, Foundation::pl_vertex);
      Foundation::MeshAttributeRegistration::execute(m, Foundation::pl_vertex);
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
      MeshControl<dim_2D>::fill_sizes(m, size_set);

      TEST_CHECK_EQUAL(size_set[0], 4u);
      TEST_CHECK_EQUAL(size_set[1], 4u);
      TEST_CHECK_EQUAL(size_set[2], 1u);

      confmeshtype_ confmesh(size_set);
      MeshControl<dim_2D>::fill_adjacencies(m, confmesh);

      typename confmeshtype_::template IndexSet<1,0>::Type& geo_vertex_at_edge(confmesh.template get_index_set<1,0>());
      typename confmeshtype_::template IndexSet<2,0>::Type& geo_vertex_at_face(confmesh.template get_index_set<2,0>());

      TEST_CHECK_EQUAL(geo_vertex_at_edge[0][0], 0u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[0][1], 1u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[1][0], 2u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[1][1], 3u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[2][0], 0u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[2][1], 2u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[3][0], 1u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[3][1], 3u);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][0], 0u);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][1], 1u);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][2], 2u);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][3], 3u);

      MeshControl<dim_2D>::fill_vertex_sets(m, confmesh, *((Attribute<double, OT_>*)(attrs.at(0).get())), *((Attribute<double, OT_>*)(attrs.at(1).get())));
      typename confmeshtype_::VertexSetType& vertex_coord_tuples(confmesh.get_vertex_set());
      TEST_CHECK_EQUAL(vertex_coord_tuples[0][0], 0u);
      TEST_CHECK_EQUAL(vertex_coord_tuples[0][1], 0u);
      TEST_CHECK_EQUAL(vertex_coord_tuples[1][0], 1u);
      TEST_CHECK_EQUAL(vertex_coord_tuples[1][1], 0u);
      TEST_CHECK_EQUAL(vertex_coord_tuples[2][0], 0u);
      TEST_CHECK_EQUAL(vertex_coord_tuples[2][1], 1u);
      TEST_CHECK_EQUAL(vertex_coord_tuples[3][0], 1u);
      TEST_CHECK_EQUAL(vertex_coord_tuples[3][1], 1u);

      delete[] size_set;
    }
};
MeshControlTest2D<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > meshcontrol_testvv("std::vector, std::vector");
MeshControlTest2D<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > meshcontrol_testvd("std::vector, std::deque");

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MeshControlTest3D:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    MeshControlTest3D(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("MeshControlTest3D<" + tag + ">")
    {
    }

    void run() const
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
      std::vector<std::shared_ptr<Foundation::AttributeBase<OT_> > > attrs;
      attrs.push_back(std::shared_ptr<Foundation::AttributeBase<OT_> >(new Foundation::Attribute<double, OT_>)); //vertex x-coords
      attrs.push_back(std::shared_ptr<Foundation::AttributeBase<OT_> >(new Foundation::Attribute<double, OT_>)); //vertex y-coords
      attrs.push_back(std::shared_ptr<Foundation::AttributeBase<OT_> >(new Foundation::Attribute<double, OT_>)); //vertex z-coords

      //8 vertices
      //face 0
      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(double(0)));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(1)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(double(0)));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(1)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(double(0)));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(1)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(1)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(double(0)));

      //complete for face 1
      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(1)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(double(1)));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(1)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(1)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(double(1)));

      //complete for face 2
      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(double(1)));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(double(0)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(double(1)));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(double(1)));

      /*  2    3
       *  *-1--*
       *  2    |
       *  |    3
       *  *--0-*
       *  0    1
       */

      //creating foundation mesh
      Foundation::Mesh<Foundation::rnt_3D, Foundation::Topology<IndexType_, OT_, IT_> > m(0, &attrs);
      Foundation::MeshAttributeRegistration::execute(m, Foundation::pl_vertex);
      Foundation::MeshAttributeRegistration::execute(m, Foundation::pl_vertex);
      Foundation::MeshAttributeRegistration::execute(m, Foundation::pl_vertex);

      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);

      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);

      m.add_polytope(Foundation::pl_face);
      m.add_polytope(Foundation::pl_face);
      m.add_polytope(Foundation::pl_face);
      m.add_polytope(Foundation::pl_face);
      m.add_polytope(Foundation::pl_face);
      m.add_polytope(Foundation::pl_face);

      m.add_polytope(Foundation::pl_polyhedron);

      //vertex 0
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 0, 0);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 0, 2);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 0, 10);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 0, 0);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 0, 3);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 0, 4);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_polyhedron, 0, 0);

      //vertex 1
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 1, 0);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 1, 3);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 1, 4);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 1, 0);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 1, 2);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 1, 4);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_polyhedron, 1, 0);

      //vertex 2
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 2, 1);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 2, 2);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 2, 11);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 2, 0);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 2, 3);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 2, 5);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_polyhedron, 2, 0);

      //vertex 3
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 3, 1);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 3, 3);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 3, 5);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 3, 2);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 3, 4);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 3, 5);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_polyhedron, 3, 0);

      //vertex 4
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 4, 4);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 4, 6);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 4, 7);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 4, 1);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 4, 2);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 4, 4);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_polyhedron, 4, 0);

      //vertex 5
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 5, 5);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 5, 6);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 5, 8);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 5, 1);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 5, 2);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 5, 5);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_polyhedron, 5, 0);

      //vertex 6
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 6, 7);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 6, 9);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 6, 10);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 6, 1);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 6, 3);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 6, 4);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_polyhedron, 6, 0);

      //vertex 7
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 7, 8);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 7, 9);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 7, 11);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 7, 1);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 7, 3);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 7, 5);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_polyhedron, 7, 0);

      typedef ConformalMesh<ConformalMeshPolicy<Shape::Hypercube<dim_3D> > > confmeshtype_;

      IndexType_* size_set(new IndexType_[4]);
      MeshControl<dim_3D>::fill_sizes(m, size_set);

      TEST_CHECK_EQUAL(size_set[0], 8u);
      TEST_CHECK_EQUAL(size_set[1], 12u);
      TEST_CHECK_EQUAL(size_set[2], 6u);
      TEST_CHECK_EQUAL(size_set[3], 1u);

      confmeshtype_ confmesh(size_set);
      MeshControl<dim_3D>::fill_adjacencies(m, confmesh);

      typename confmeshtype_::template IndexSet<1,0>::Type& geo_vertex_at_edge(confmesh.template get_index_set<1,0>());
      typename confmeshtype_::template IndexSet<2,0>::Type& geo_vertex_at_face(confmesh.template get_index_set<2,0>());
      typename confmeshtype_::template IndexSet<3,0>::Type& geo_vertex_at_poly(confmesh.template get_index_set<3,0>());

      //TEST_CHECK_EQUAL(geo_vertex_at_edge[0][0], 0u);

      MeshControl<dim_3D>::fill_vertex_sets(m, confmesh, *((Attribute<double, OT_>*)(attrs.at(0).get())), *((Attribute<double, OT_>*)(attrs.at(1).get())), *((Attribute<double, OT_>*)(attrs.at(2).get())));

      typename confmeshtype_::VertexSetType& vertex_coord_tuples(confmesh.get_vertex_set());
      //TEST_CHECK_EQUAL(vertex_coord_tuples[0][0], 0u);

      delete[] size_set;
    }
};
MeshControlTest3D<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > meshcontrol3d_testvv("std::vector, std::vector");
MeshControlTest3D<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > meshcontrol3d_testvd("std::vector, std::deque");
