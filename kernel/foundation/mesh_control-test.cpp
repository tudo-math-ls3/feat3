#define SERIAL
#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/mesh_control.hpp>
#include<kernel/foundation/dense_data_wrapper.hpp>
#include<kernel/geometry/macro_factory.hpp>
#include<kernel/lafem/dense_vector.hpp>
#include<kernel/lafem/vector_mirror.hpp>
#include<kernel/foundation/halo_control.hpp>
#include<kernel/foundation/halo.hpp>
#include<kernel/geometry/cell_sub_set.hpp>
#include<kernel/archs.hpp>
#include<deque>
#include<algorithm>

#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/dof_adjacency.hpp>
#include <kernel/space/dof_mirror.hpp>
#include <kernel/assembly/standard_operators.hpp>
#include <kernel/assembly/standard_functionals.hpp>
#include <kernel/assembly/dirichlet_bc.hpp>

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

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(1));

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

      typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

      IndexType_* size_set(new IndexType_[2]);
      MeshControl<dim_1D>::fill_sizes(m, size_set);

      TEST_CHECK_EQUAL(size_set[0], Index(2));
      TEST_CHECK_EQUAL(size_set[1], Index(1));

      confmeshtype_ confmesh(size_set);
      MeshControl<dim_1D>::fill_adjacencies(m, confmesh);

      typename confmeshtype_::template IndexSet<1,0>::Type& geo_vertex_at_edge(confmesh.template get_index_set<1,0>());
      TEST_CHECK_EQUAL(geo_vertex_at_edge[0][0], 0u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[0][1], 1u);

      MeshControl<dim_1D>::fill_vertex_sets(m, confmesh, *((Attribute<double, OT_>*)(attrs.at(0).get())));
      typename confmeshtype_::VertexSetType& vertex_coord_tuples(confmesh.get_vertex_set());
      TEST_CHECK_EQUAL(vertex_coord_tuples[0][0], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[1][0], double(1));

      ///test reconversion
      Foundation::Mesh<Foundation::rnt_1D, Foundation::Topology<IndexType_, OT_, IT_> > m1(1);
      MeshControl<dim_1D>::fill_adjacencies(confmesh, m1);

      for(Index i(0) ; i < m.get_topologies().size() ; ++i)
        for(Index j(0) ; j < m.get_topologies().at(i).size() ; ++j)
          for(Index k(0) ; k < m.get_topologies().at(i).at(j).size() ; ++k)
            TEST_CHECK_EQUAL(m1.get_topologies().at(i).at(j).at(k), m.get_topologies().at(i).at(j).at(k));

      MeshControl<dim_1D>::fill_vertex_sets(confmesh, m1, *((Attribute<double, OT_>*)(attrs.at(0).get())));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(0), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(1), double(1));

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

      typedef ConformalMesh<Shape::Hypercube<dim_2D> > confmeshtype_;

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

      ///test reconversion
      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m1(1);
      MeshControl<dim_2D>::fill_adjacencies(confmesh, m1);

      for(Index i(0) ; i < m.get_topologies().size() ; ++i)
        for(Index j(0) ; j < m.get_topologies().at(i).size() ; ++j)
          for(Index k(0) ; k < m.get_topologies().at(i).at(j).size() ; ++k)
            TEST_CHECK_EQUAL(m1.get_topologies().at(i).at(j).at(k), m.get_topologies().at(i).at(j).at(k));

      MeshControl<dim_2D>::fill_vertex_sets(confmesh, m1, *((Attribute<double, OT_>*)(attrs.at(0).get())), *((Attribute<double, OT_>*)(attrs.at(1).get())));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(0), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(0), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(1), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(1), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(2), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(2), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(3), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(3), double(1));


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
      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(0));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(1));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(0));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(1));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(0));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(1));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(1));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(0));

      //complete for face 1
      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(1));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(1));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(1));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(1));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(1));

      //complete for face 2
      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(1));

      ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().push_back(double(0));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().push_back(double(1));
      ((Foundation::Attribute<double, OT_>*)(attrs.at(2).get()))->get_data().push_back(double(1));

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
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_face, 3, 0);
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

      typedef ConformalMesh<Shape::Hypercube<dim_3D> > confmeshtype_;

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

      TEST_CHECK_EQUAL(geo_vertex_at_edge[0][0], 0u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[0][1], 1u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[1][0], 2u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[1][1], 3u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[2][0], 0u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[2][1], 2u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[3][0], 1u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[3][1], 3u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[4][0], 1u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[4][1], 4u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[5][0], 3u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[5][1], 5u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[6][0], 4u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[6][1], 5u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[7][0], 4u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[7][1], 6u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[8][0], 5u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[8][1], 7u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[9][0], 6u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[9][1], 7u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[10][0], 0u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[10][1], 6u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[11][0], 2u);
      TEST_CHECK_EQUAL(geo_vertex_at_edge[11][1], 7u);

      TEST_CHECK_EQUAL(geo_vertex_at_face[0][0], 0ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][1], 1ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][2], 2ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[0][3], 3ul);

      TEST_CHECK_EQUAL(geo_vertex_at_face[1][0], 4ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[1][1], 5ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[1][2], 6ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[1][3], 7ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[2][0], 1ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[2][1], 3ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[2][2], 4ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[2][3], 5ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[3][0], 0ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[3][1], 2ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[3][2], 6ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[3][3], 7ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[4][0], 0ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[4][1], 1ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[4][2], 4ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[4][3], 6ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[5][0], 2ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[5][1], 3ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[5][2], 5ul);
      TEST_CHECK_EQUAL(geo_vertex_at_face[5][3], 7ul);

      TEST_CHECK_EQUAL(geo_vertex_at_poly[0][0], 0ul);
      TEST_CHECK_EQUAL(geo_vertex_at_poly[0][1], 1ul);
      TEST_CHECK_EQUAL(geo_vertex_at_poly[0][2], 2ul);
      TEST_CHECK_EQUAL(geo_vertex_at_poly[0][3], 3ul);
      TEST_CHECK_EQUAL(geo_vertex_at_poly[0][4], 4ul);
      TEST_CHECK_EQUAL(geo_vertex_at_poly[0][5], 5ul);
      TEST_CHECK_EQUAL(geo_vertex_at_poly[0][6], 6ul);
      TEST_CHECK_EQUAL(geo_vertex_at_poly[0][7], 7ul);

      MeshControl<dim_3D>::fill_vertex_sets(m, confmesh, *((Attribute<double, OT_>*)(attrs.at(0).get())), *((Attribute<double, OT_>*)(attrs.at(1).get())), *((Attribute<double, OT_>*)(attrs.at(2).get())));

      typename confmeshtype_::VertexSetType& vertex_coord_tuples(confmesh.get_vertex_set());
      TEST_CHECK_EQUAL(vertex_coord_tuples[0][0], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[0][1], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[0][2], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[1][0], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[1][1], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[1][2], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[2][0], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[2][1], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[2][2], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[3][0], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[3][1], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[3][2], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[4][0], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[4][1], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[4][2], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[5][0], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[5][1], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[5][2], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[6][0], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[6][1], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[6][2], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[7][0], double(0));
      TEST_CHECK_EQUAL(vertex_coord_tuples[7][1], double(1));
      TEST_CHECK_EQUAL(vertex_coord_tuples[7][2], double(1));

      ///test reconversion
      Foundation::Mesh<Foundation::rnt_3D, Foundation::Topology<IndexType_, OT_, IT_> > m1(1);
      MeshControl<dim_3D>::fill_adjacencies(confmesh, m1);

      MeshControl<dim_3D>::fill_vertex_sets(confmesh, m1, *((Attribute<double, OT_>*)(attrs.at(0).get())), *((Attribute<double, OT_>*)(attrs.at(1).get())), *((Attribute<double, OT_>*)(attrs.at(2).get())));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(0), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(0), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(2).get()))->at(0), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(1), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(1), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(2).get()))->at(1), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(2), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(2), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(2).get()))->at(2), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(3), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(3), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(2).get()))->at(3), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(4), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(4), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(2).get()))->at(4), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(5), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(5), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(2).get()))->at(5), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(6), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(6), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(2).get()))->at(6), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(0).get()))->at(7), double(0));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(1).get()))->at(7), double(1));
      TEST_CHECK_EQUAL(((Attribute<double, OT_>*)(attrs.at(2).get()))->at(7), double(1));

      delete[] size_set;
    }
};
MeshControlTest3D<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > meshcontrol3d_testvv("std::vector, std::vector");
MeshControlTest3D<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > meshcontrol3d_testvd("std::vector, std::deque");

template<typename T_>
class RhsFunc
{
public:
  static T_ eval(T_ /*x*/, T_ /*y*/)
  {
    return T_(1);
  }
};
template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MeshControlPartitioningTest2D:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    MeshControlPartitioningTest2D(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("MeshControlPartitioningTest2D<" + tag + ">")
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

      ///get a conformal mesh as basemesh
      typedef ConformalMesh<Shape::Hypercube<dim_2D> > basemeshtype_;

      IndexType_* size_set(new IndexType_[3]);
      MeshControl<dim_2D>::fill_sizes(m, size_set);

      basemeshtype_ basemesh(size_set);
      MeshControl<dim_2D>::fill_adjacencies(m, basemesh);
      MeshControl<dim_2D>::fill_vertex_sets(m, basemesh, *((Attribute<double, OT_>*)(attrs.at(0).get())), *((Attribute<double, OT_>*)(attrs.at(1).get())));

      ///refine basemesh to match process count
      Geometry::StandardRefinery<basemeshtype_> mesh_refinery(basemesh);
      basemeshtype_ fine_basemesh(mesh_refinery);

      ///get foundation mesh, create halo set
      Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > fine_bm_found(4711, &attrs);
      MeshControl<dim_2D>::fill_adjacencies(fine_basemesh, fine_bm_found);
      MeshControl<dim_2D>::fill_vertex_sets(fine_basemesh, fine_bm_found, *((Attribute<double, OT_>*)(attrs.at(0).get())), *((Attribute<double, OT_>*)(attrs.at(1).get())));

      ///depending on rank: compute adjacent macros to potentially communicate with, e.g. rank 0
      typedef typename Topology<IndexType_, OT_, IT_>::storage_type_ TopologyStorageType;
      TopologyStorageType potential_comm_partners_for_face_0(fine_bm_found.get_adjacent_polytopes(pl_face, pl_face, 0));

      ///create halos
      std::vector<std::shared_ptr<HaloBase<Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > > > > halos;
      for(Index i(0) ; i < potential_comm_partners_for_face_0.size() ; ++i)
      {
        TopologyStorageType comm_intersect_0_i(fine_bm_found.get_comm_intersection(pl_face, pl_edge, 0, i));
        if(comm_intersect_0_i.size() == 0)
        {
          comm_intersect_0_i = fine_bm_found.get_comm_intersection(pl_face, pl_vertex, 0, i);
          for(Index j(0) ; j < comm_intersect_0_i.size() ; ++j)
          {
            Halo<0, pl_vertex, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > > halo(fine_bm_found, i);
            halo.add_element_pair(comm_intersect_0_i.at(j), comm_intersect_0_i.at(j));

            halos.push_back(std::shared_ptr<HaloBase<Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > > >(new Halo<0, pl_vertex, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > >(halo)));
          }
        }
        else
        {
          for(Index j(0) ; j < comm_intersect_0_i.size() ; ++j)
          {
            Halo<0, pl_edge, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > > halo(fine_bm_found, i);
            halo.add_element_pair(comm_intersect_0_i.at(j), comm_intersect_0_i.at(j));
            halos.push_back(std::shared_ptr<HaloBase<Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > > >(new Halo<0, pl_edge, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > >(halo)));
          }
        }
      }
      std::sort(halos.begin(), halos.end(), compare_other<Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> >, OT_, IndexType_>);

      TEST_CHECK_EQUAL(halos.size(), 3);
      TEST_CHECK_EQUAL(halos.at(0)->get_other(), 1);
      TEST_CHECK_EQUAL(halos.at(1)->get_other(), 2);
      TEST_CHECK_EQUAL(halos.at(2)->get_other(), 3);

      ///before refinement, get the actual macro for assembly (of inner BC mainly)
      Index macro_number(0);
      MacroFactory<basemeshtype_> macro_factory(fine_basemesh, macro_number);
      basemeshtype_ macro_mesh(macro_factory);
      Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > macro_mesh_found(99);
      MeshControl<dim_2D>::fill_adjacencies(macro_mesh, macro_mesh_found);
      Attribute<double, OT_> vertex_x_coords_macro;
      Attribute<double, OT_> vertex_y_coords_macro;
      MeshControl<dim_2D>::fill_vertex_sets(macro_mesh, macro_mesh_found, vertex_x_coords_macro, vertex_y_coords_macro);

      TEST_CHECK_EQUAL(vertex_x_coords_macro.at(0), double(0));
      TEST_CHECK_EQUAL(vertex_y_coords_macro.at(0), double(0));
      TEST_CHECK_EQUAL(vertex_x_coords_macro.at(1), double(0.5));
      TEST_CHECK_EQUAL(vertex_y_coords_macro.at(1), double(0));
      TEST_CHECK_EQUAL(vertex_x_coords_macro.at(2), double(0));
      TEST_CHECK_EQUAL(vertex_y_coords_macro.at(2), double(0.5));
      TEST_CHECK_EQUAL(vertex_x_coords_macro.at(3), double(0.5));
      TEST_CHECK_EQUAL(vertex_y_coords_macro.at(3), double(0.5));

      //every edge is an inner boundary
      std::vector<std::shared_ptr<HaloBase<Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > > > > macro_boundaries;

      for(Index i(0) ; i < macro_mesh_found.get_topologies().at(ipi_edge_vertex).size() ; ++i)
      {
        Halo<0, pl_edge, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > > result(macro_mesh_found);
        result.add_element_pair(i, i);
        macro_boundaries.push_back(std::shared_ptr<HaloBase<Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > > >(new Halo<0, pl_edge, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > >(result)));
      }
      TEST_CHECK_EQUAL(macro_boundaries.size(), 4);

      ///refine basemesh, macromesh
      Geometry::StandardRefinery<basemeshtype_> mesh_refinery_fine(fine_basemesh);
      Geometry::StandardRefinery<basemeshtype_> mesh_refinery_macro(macro_mesh);
      basemeshtype_ finemost_basemesh(mesh_refinery_fine);
      basemeshtype_ finemost_macro(mesh_refinery_macro);

      ///get refined halos
      std::vector<std::shared_ptr<Geometry::CellSubSet<Shape::Hypercube<2> > > > finemost_cell_subsets;
      for(Index i(0) ; i < halos.size() ; ++i)
      {
        //transform found->geo
        Index* polytopes_in_subset = new Index[3];

        if(halos.at(i)->get_level() == pl_vertex)
          HaloControl<dim_1D>::fill_sizes(*((Halo<0, pl_vertex, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > >*)(halos.at(i).get())), polytopes_in_subset);
        else
          HaloControl<dim_2D>::fill_sizes(*((Halo<0, pl_edge, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > >*)(halos.at(i).get())), polytopes_in_subset);

        Geometry::CellSubSet<Shape::Hypercube<2> > cell_sub_set(polytopes_in_subset);

        if(halos.at(i)->get_level() == pl_vertex)
          HaloControl<dim_1D>::fill_target_set(*((Halo<0, pl_vertex, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > >*)(halos.at(i).get())), cell_sub_set);
        else
          HaloControl<dim_2D>::fill_target_set(*((Halo<0, pl_edge, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > >*)(halos.at(i).get())), cell_sub_set);

        //refine
        Geometry::StandardRefinery<Geometry::CellSubSet<Shape::Hypercube<2> >, basemeshtype_> cell_refinery(cell_sub_set, fine_basemesh);

        //add
        finemost_cell_subsets.push_back(std::shared_ptr<Geometry::CellSubSet<Shape::Hypercube<2> > >(new Geometry::CellSubSet<Shape::Hypercube<2> >(cell_refinery)));

        delete[] polytopes_in_subset;
      }

      ///get refined macro boundaries
      std::vector<std::shared_ptr<Geometry::CellSubSet<Shape::Hypercube<2> > > > finemost_macro_boundaries;
      for(Index i(0) ; i < macro_boundaries.size() ; ++i)
      {
        //transform found->geo
        Index* polytopes_in_subset = new Index[3];
        HaloControl<dim_2D>::fill_sizes(*((Halo<0, pl_edge, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > >*)(macro_boundaries.at(i).get())), polytopes_in_subset);
        Geometry::CellSubSet<Shape::Hypercube<2> > cell_sub_set(polytopes_in_subset);
        HaloControl<dim_2D>::fill_target_set(*((Halo<0, pl_edge, Mesh<rnt_2D, Topology<IndexType_, OT_, IT_> > >*)(macro_boundaries.at(i).get())), cell_sub_set);

        //refine
        Geometry::StandardRefinery<Geometry::CellSubSet<Shape::Hypercube<2> >, basemeshtype_> cell_refinery(cell_sub_set, macro_mesh);

        //add
        finemost_macro_boundaries.push_back(std::shared_ptr<Geometry::CellSubSet<Shape::Hypercube<2> > >(new Geometry::CellSubSet<Shape::Hypercube<2> >(cell_refinery)));

        delete[] polytopes_in_subset;
      }
      TEST_CHECK_EQUAL(finemost_macro_boundaries.size(), 4);

      ///assembly
      ///we now have: finemost basemesh, macro_mesh and its boundary components, cell_subsets from halos all in 'geometry mode'
      // create trafo
      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > trafo(finemost_macro);
      // create space
      Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);

      SparseMatrixCSR<Mem::Main, double> mat_sys(Space::DofAdjacency<>::assemble(space));
      mat_sys.clear();
      Assembly::BilinearScalarLaplaceFunctor::assemble(mat_sys, space, "gauss-legendre:2");

      DenseVector<Mem::Main, double> vec_rhs(space.get_num_dofs(), double(0));
      Assembly::LinearScalarIntegralFunctor<RhsFunc>::assemble(vec_rhs, space, "gauss-legendre:2");

      // assemble homogeneous Dirichlet BCs
      Assembly::DirichletBC<Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > > dirichlet(space);
      for(Index i(0) ; i < finemost_macro_boundaries.size() ; ++i)
      {
        dirichlet.add_cell_set(*finemost_macro_boundaries.at(i).get());
      }

      // allocate solution vector
      DenseVector<Mem::Main, double> vec_sol(space.get_num_dofs(), double(0));

      // assemble filter:
      UnitFilter<Mem::Main, double> filter(dirichlet.template assemble<Mem::Main, double>());

      // filter system
      filter.filter_mat(mat_sys);
      filter.filter_rhs(vec_rhs);
      filter.filter_sol(vec_sol);

      ///assemble mirrors
      std::vector<LAFEM::VectorMirror<Mem::Main, double> > mirrors;
      std::vector<LAFEM::DenseVector<Mem::Main, double> > sendbufs;
      std::vector<LAFEM::DenseVector<Mem::Main, double> > recvbufs;
      std::vector<Index> destranks;
      std::vector<Index> sourceranks;

      Graph dof_adj(Space::DofAdjacency<>::assemble(space));

      for(Index i(0) ; i < finemost_macro_boundaries.size() ; ++i)
      {
        Graph dof_mirror(Space::DofMirror::assemble(space, *(finemost_macro_boundaries.at(i).get())));
        VectorMirror<Mem::Main, double> target_mirror(dof_mirror);
        DenseVector<Mem::Main, double> sendbuf(target_mirror.size());
        DenseVector<Mem::Main, double> recvbuf(target_mirror.size());

        mirrors.push_back(target_mirror);
        sendbufs.push_back(sendbuf);
        recvbufs.push_back(recvbuf);

        destranks.push_back(macro_boundaries.at(i)->get_other());
        sourceranks.push_back(macro_number);
      }
      delete[] size_set;

    }
};
MeshControlPartitioningTest2D<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > meshcontrolpart_testvv("std::vector, std::vector");
MeshControlPartitioningTest2D<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > meshcontrolpart_testvd("std::vector, std::deque");
