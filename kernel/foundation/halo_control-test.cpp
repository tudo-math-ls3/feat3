#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/halo_control.hpp>
#include<kernel/foundation/dense_data_wrapper.hpp>
#include<kernel/geometry/cell_sub_set.hpp>
#include<kernel/lafem/dense_vector.hpp>
#include<kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;
using namespace FEAST::Geometry;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class HaloControlTest1D:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    HaloControlTest1D(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("HaloControlTest1D<" + tag + ">")
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

      //creating conformal geometry mesh
      IndexType_ index_set[] = {2, 1};
      typedef Geometry::ConformalMesh<Geometry::ConformalMeshPolicy<Shape::Hypercube<1> > > confmeshtype_;
      confmeshtype_ geo_m(index_set);

      //transfer data
      //Edge->Vertex
      typename confmeshtype_::template IndexSet<1,0>::Type& geo_vertex_at_edge(geo_m.template get_index_set<1,0>());
      for(IndexType_ i(0) ; i < m.get_topologies().at(Foundation::ipi_edge_vertex).size() ; ++i)
      {
        //get all adjacencies Edge->Vertex from foundation mesh
        typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ found_vertex_at_edge_i(m.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, i));

        for(IndexType_ j(0) ; j < found_vertex_at_edge_i.size() ; ++j)
        {
          geo_vertex_at_edge[i][j] = found_vertex_at_edge_i.at(j); //edge i, adjacent vertex j
        }
      }

      typename confmeshtype_::VertexSetType& vertex_coord_tuples(geo_m.get_vertex_set());
      vertex_coord_tuples[0][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(0); //xcoord of first node
      vertex_coord_tuples[1][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(1);

      //-------------------------------------
      //create halo with one vertex-vertex pair
      Foundation::Halo<0, Foundation::pl_vertex, Foundation::Mesh<Foundation::rnt_1D, Foundation::Topology<IndexType_, OT_, IT_> > > h(m, 1);
      h.add_element_pair(1, 0);

      //create CellSubSet
      Index* polytopes_in_subset = new Index[2]; //no overlap and one vertex means one vertex
      Foundation::HaloControl<Foundation::dim_1D>::fill_sizes(h, polytopes_in_subset);

      TEST_CHECK_EQUAL(polytopes_in_subset[0], Index(1));
      TEST_CHECK_EQUAL(polytopes_in_subset[1], Index(0));

      CellSubSet<Shape::Hypercube<1> > cell_sub_set(polytopes_in_subset);
      HaloControl<Foundation::dim_1D>::fill_target_set(h, cell_sub_set);

      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[0], 1ul);
      //-------------------------------------

      //create halo with one edge-edge pair
      Foundation::Halo<1, Foundation::pl_edge, Foundation::Mesh<Foundation::rnt_1D, Foundation::Topology<IndexType_, OT_, IT_> > > h1(m, 1);
      h1.add_element_pair(0, 0);

      Index* polytopes_in_subset1 = new Index[2]; //overlap 1 and one edge means 2 vertices
      Foundation::HaloControl<Foundation::dim_1D>::fill_sizes(h1, polytopes_in_subset1);

      TEST_CHECK_EQUAL(polytopes_in_subset1[0], Index(2));
      TEST_CHECK_EQUAL(polytopes_in_subset1[1], Index(1));

      CellSubSet<Shape::Hypercube<1> > cell_sub_set1(polytopes_in_subset1);
      HaloControl<Foundation::dim_1D>::fill_target_set(h1, cell_sub_set1);

      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[0], 0ul); //edge 0
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[0], 0ul); //vertex 0
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[1], 1ul); //vertex 1
      //-------------------------------------

      delete[] polytopes_in_subset;
      delete[] polytopes_in_subset1;
    }
};
HaloControlTest1D<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > halocontrol1d_testvv("std::vector, std::vector");
HaloControlTest1D<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > halocontrol1d_testvd("std::vector, std::deque");



template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class HaloControlTest2D:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    HaloControlTest2D(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("HaloControlTest2D<" + tag + ">")
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


      //creating conformal geometry mesh
      IndexType_ index_set[] = {4, 4, 1};
      typedef Geometry::ConformalMesh<Geometry::ConformalMeshPolicy<Shape::Hypercube<2> > > confmeshtype_;
      confmeshtype_ geo_m(index_set);

      //transfer data
      //Edge->Vertex
      typename confmeshtype_::template IndexSet<1,0>::Type& geo_vertex_at_edge(geo_m.template get_index_set<1,0>());
      for(IndexType_ i(0) ; i < m.get_topologies().at(Foundation::ipi_edge_vertex).size() ; ++i)
      {
        //get all adjacencies Edge->Vertex from foundation mesh
        typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ found_vertex_at_edge_i(m.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, i));

        for(IndexType_ j(0) ; j < found_vertex_at_edge_i.size() ; ++j)
        {
          geo_vertex_at_edge[i][j] = found_vertex_at_edge_i.at(j); //edge i, adjacent vertex j
        }
      }
      //Face->Vertex
      typename confmeshtype_::template IndexSet<2,0>::Type& geo_vertex_at_face(geo_m.template get_index_set<2,0>());
      for(IndexType_ i(0) ; i < m.get_topologies().at(Foundation::ipi_face_vertex).size() ; ++i)
      {
        //get all adjacencies Face->Vertex from foundation mesh
        typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ found_vertex_at_face_i(m.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, i));

        for(IndexType_ j(0) ; j < found_vertex_at_face_i.size() ; ++j)
        {
          geo_vertex_at_face[i][j] = found_vertex_at_face_i.at(j); //face i, adjacent vertex j
        }
      }

      typename confmeshtype_::VertexSetType& vertex_coord_tuples(geo_m.get_vertex_set());
      vertex_coord_tuples[0][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(0); //xcoord of first node
      vertex_coord_tuples[0][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(0); //ycoord of first node
      vertex_coord_tuples[1][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(1);
      vertex_coord_tuples[1][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(1);
      vertex_coord_tuples[2][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(2);
      vertex_coord_tuples[2][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(2);
      vertex_coord_tuples[3][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(3);
      vertex_coord_tuples[3][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(3);

      //-----------------------------------------
      //create halo with one edge-edge pair
      Foundation::Halo<0, Foundation::pl_edge, Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > > h(m, 1);
      h.add_element_pair(3, 2);

      //create CellSubSet
      Index* polytopes_in_subset = new Index[3]; //no overlap and one edge means two vertices but no faces
      Foundation::HaloControl<Foundation::dim_2D>::fill_sizes(h, polytopes_in_subset);

      TEST_CHECK_EQUAL(polytopes_in_subset[0], Index(2));
      TEST_CHECK_EQUAL(polytopes_in_subset[1], Index(1));
      TEST_CHECK_EQUAL(polytopes_in_subset[2], Index(0));

      Geometry::CellSubSet<Shape::Hypercube<2> > cell_sub_set(polytopes_in_subset);
      Foundation::HaloControl<Foundation::dim_2D>::fill_target_set(h, cell_sub_set);

      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[0], 1ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[1], 3ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<1>()[0], 3ul);

      //-----------------------------------------
      //create halo with one face-face pair
      Foundation::Halo<1, Foundation::pl_face, Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > > h1(m, 1);
      h1.add_element_pair(0, 0);

      Index* polytopes_in_subset1 = new Index[3]; //overlap 1 and one face means 4 vertices and 4 edges
      Foundation::HaloControl<Foundation::dim_2D>::fill_sizes(h1, polytopes_in_subset1);

      TEST_CHECK_EQUAL(polytopes_in_subset1[0], Index(4));
      TEST_CHECK_EQUAL(polytopes_in_subset1[1], Index(4));
      TEST_CHECK_EQUAL(polytopes_in_subset1[2], Index(1));

      Geometry::CellSubSet<Shape::Hypercube<2> > cell_sub_set1(polytopes_in_subset1);
      Foundation::HaloControl<Foundation::dim_2D>::fill_target_set(h1, cell_sub_set1);

      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[0], 0ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[1], 1ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[2], 2ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[3], 3ul);

      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[0], 0ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[1], 2ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[2], 3ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[3], 1ul);

      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<2>()[0], 0ul);

      //-----------------------------------------

      delete[] polytopes_in_subset;
      delete[] polytopes_in_subset1;
    }
};
HaloControlTest2D<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > halocontrol_testvv("std::vector, std::vector");
HaloControlTest2D<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > halocontrol_testvd("std::vector, std::deque");

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class HaloControlTest3D:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    HaloControlTest3D(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("HaloControlTest3D<" + tag + ">")
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

      //creating conformal geometry mesh
      IndexType_ index_set[] = {8, 12, 6, 1};
      typedef Geometry::ConformalMesh<Geometry::ConformalMeshPolicy<Shape::Hypercube<3> > > confmeshtype_;
      confmeshtype_ geo_m(index_set);

      //transfer data
      //Edge->Vertex
      typename confmeshtype_::template IndexSet<1,0>::Type& geo_vertex_at_edge(geo_m.template get_index_set<1,0>());
      for(IndexType_ i(0) ; i < m.get_topologies().at(Foundation::ipi_edge_vertex).size() ; ++i)
      {
        //get all adjacencies Edge->Vertex from foundation mesh
        typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ found_vertex_at_edge_i(m.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, i));

        for(IndexType_ j(0) ; j < found_vertex_at_edge_i.size() ; ++j)
        {
          geo_vertex_at_edge[i][j] = found_vertex_at_edge_i.at(j); //edge i, adjacent vertex j
        }
      }
      //Face->Vertex
      typename confmeshtype_::template IndexSet<2,0>::Type& geo_vertex_at_face(geo_m.template get_index_set<2,0>());
      for(IndexType_ i(0) ; i < m.get_topologies().at(Foundation::ipi_face_vertex).size() ; ++i)
      {
        //get all adjacencies Face->Vertex from foundation mesh
        typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ found_vertex_at_face_i(m.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, i));

        for(IndexType_ j(0) ; j < found_vertex_at_face_i.size() ; ++j)
        {
          geo_vertex_at_face[i][j] = found_vertex_at_face_i.at(j); //face i, adjacent vertex j
        }
      }
      //Poly->Vertex
      typename confmeshtype_::template IndexSet<3,0>::Type& geo_vertex_at_polyhedron(geo_m.template get_index_set<3,0>());
      for(IndexType_ i(0) ; i < m.get_topologies().at(Foundation::ipi_polyhedron_vertex).size() ; ++i)
      {
        //get all adjacencies poly->vertex from foundation mesh
        typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ found_vertex_at_polyhedron_i(m.get_adjacent_polytopes(Foundation::pl_polyhedron, Foundation::pl_vertex, i));

        for(IndexType_ j(0) ; j < found_vertex_at_polyhedron_i.size() ; ++j)
        {
          geo_vertex_at_polyhedron[i][j] = found_vertex_at_polyhedron_i.at(j); //poly i, adjacent vertex j
        }
      }

      typename confmeshtype_::VertexSetType& vertex_coord_tuples(geo_m.get_vertex_set());
      vertex_coord_tuples[0][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(0); //xcoord of first node
      vertex_coord_tuples[0][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(0); //ycoord of first node
      vertex_coord_tuples[0][2] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(2).get()))->get_data().at(0); //zcoord of first node
      vertex_coord_tuples[1][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(1);
      vertex_coord_tuples[1][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(1);
      vertex_coord_tuples[1][2] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(2).get()))->get_data().at(1);
      vertex_coord_tuples[2][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(2);
      vertex_coord_tuples[2][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(2);
      vertex_coord_tuples[2][2] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(2).get()))->get_data().at(2);
      vertex_coord_tuples[3][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(3);
      vertex_coord_tuples[3][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(3);
      vertex_coord_tuples[3][2] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(2).get()))->get_data().at(3);
      vertex_coord_tuples[4][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(4);
      vertex_coord_tuples[4][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(4);
      vertex_coord_tuples[4][2] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(2).get()))->get_data().at(4);
      vertex_coord_tuples[5][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(5);
      vertex_coord_tuples[5][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(5);
      vertex_coord_tuples[5][2] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(2).get()))->get_data().at(5);
      vertex_coord_tuples[6][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(6);
      vertex_coord_tuples[6][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(6);
      vertex_coord_tuples[6][2] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(2).get()))->get_data().at(6);
      vertex_coord_tuples[7][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(7);
      vertex_coord_tuples[7][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(7);
      vertex_coord_tuples[7][2] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(2).get()))->get_data().at(7);


      //create halo with one face-face pair
      Foundation::Halo<0, Foundation::pl_face, Foundation::Mesh<Foundation::rnt_3D, Foundation::Topology<IndexType_, OT_, IT_> > > h(m, 1);
      h.add_element_pair(2, 3);

      //create CellSubSet
      Index* polytopes_in_subset = new Index[4]; //no overlap and one face means four vertices and four edges
      Foundation::HaloControl<Foundation::dim_3D>::fill_sizes(h, polytopes_in_subset);

      TEST_CHECK_EQUAL(polytopes_in_subset[0], Index(4));
      TEST_CHECK_EQUAL(polytopes_in_subset[1], Index(4));
      TEST_CHECK_EQUAL(polytopes_in_subset[2], Index(1));
      TEST_CHECK_EQUAL(polytopes_in_subset[3], Index(0));

      CellSubSet<Shape::Hypercube<3> > cell_sub_set(polytopes_in_subset);
      Foundation::HaloControl<Foundation::dim_3D>::fill_target_set(h, cell_sub_set);

      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[0], 1ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[1], 3ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[2], 4ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[3], 5ul);

      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<1>()[0], 3ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<1>()[1], 4ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<1>()[2], 5ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<1>()[3], 6ul);

      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<2>()[0], 2ul);

      //create halo with one poly-poly pair
      Foundation::Halo<1, Foundation::pl_polyhedron, Foundation::Mesh<Foundation::rnt_3D, Foundation::Topology<IndexType_, OT_, IT_> > > h1(m, 1);
      h1.add_element_pair(0, 0);

      Index* polytopes_in_subset1 = new Index[4]; //overlap 1 and one face means 4 vertices and 4 edges
      Foundation::HaloControl<Foundation::dim_3D>::fill_sizes(h1, polytopes_in_subset1);

      TEST_CHECK_EQUAL(polytopes_in_subset1[0], Index(8));
      TEST_CHECK_EQUAL(polytopes_in_subset1[1], Index(12));
      TEST_CHECK_EQUAL(polytopes_in_subset1[2], Index(6));
      TEST_CHECK_EQUAL(polytopes_in_subset1[3], Index(1));

      CellSubSet<Shape::Hypercube<3> > cell_sub_set1(polytopes_in_subset1);
      Foundation::HaloControl<Foundation::dim_3D>::fill_target_set(h1, cell_sub_set1);

      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[0], 0ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[1], 1ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[2], 2ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[3], 6ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[4], 7ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[5], 3ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[6], 4ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<0>()[7], 5ul);

      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[0], 0ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[1], 2ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[2], 10ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[3], 11ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[4], 9ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[5], 3ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[6], 4ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[7], 7ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[8], 5ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[9], 6ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[10], 1ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<1>()[11], 8ul);

      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<2>()[0], 0ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<2>()[1], 3ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<2>()[2], 4ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<2>()[3], 2ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<2>()[4], 5ul);
      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<2>()[5], 1ul);

      TEST_CHECK_EQUAL(cell_sub_set1.template get_target_set<3>()[0], 0ul);

      delete[] polytopes_in_subset;
      delete[] polytopes_in_subset1;
    }
};
HaloControlTest3D<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > halocontrol3d_testvv("std::vector, std::vector");
HaloControlTest3D<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > halocontrol3d_testvd("std::vector, std::deque");
