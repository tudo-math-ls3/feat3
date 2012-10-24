#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/halo_control.hpp>
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
class HaloControlTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    HaloControlTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("HaloControlTest<" + tag + ">")
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

      //create halo with one edge-edge pair
      Foundation::Halo<0, Foundation::pl_edge, Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > > h(m, 1);
      h.add_element_pair(3, 2);

      //create CellSubSet
      Index* polytopes_in_subset = new Index[3]; //no overlap and one edge means two vertices but no faces
      Foundation::HaloControl<Foundation::dim_2D>::fill_sizes(h, polytopes_in_subset);

      TEST_CHECK_EQUAL(polytopes_in_subset[0], Index(2));
      TEST_CHECK_EQUAL(polytopes_in_subset[1], Index(1));
      TEST_CHECK_EQUAL(polytopes_in_subset[2], Index(0));

      //create halo with one face-face pair
      Foundation::Halo<1, Foundation::pl_face, Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > > h1(m, 1);
      h1.add_element_pair(0, 0);

      Index* polytopes_in_subset1 = new Index[3]; //no overlap and one edge means two vertices but no faces
      Foundation::HaloControl<Foundation::dim_2D>::fill_sizes(h1, polytopes_in_subset1);

      TEST_CHECK_EQUAL(polytopes_in_subset1[0], Index(4));
      TEST_CHECK_EQUAL(polytopes_in_subset1[1], Index(4));
      TEST_CHECK_EQUAL(polytopes_in_subset1[2], Index(1));
      /*Geometry::CellSubSet<Shape::Hypercube<2> > cell_sub_set(polytopes_in_subset);

      typename Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> >::storage_type_ adjacent_vertices(h.get_mesh().get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, h.get_element(0)));

      cell_sub_set.template get_target_set<0>()[0] = adjacent_vertices.at(0); //first vertex
      cell_sub_set.template get_target_set<0>()[1] = adjacent_vertices.at(1); //second vertex
      cell_sub_set.template get_target_set<1>()[0] = h.get_element(0); //only edge

      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[0], 1ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[1], 3ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<1>()[0], 3ul);*/

      delete[] polytopes_in_subset;
    }
};
HaloControlTest<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > halocontrol_testvv("std::vector, std::vector");
HaloControlTest<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > halocontrol_testvd("std::vector, std::deque");
