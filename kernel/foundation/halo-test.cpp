#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/base.hpp>
#include <kernel/foundation/attribute.hpp>
#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/archs.hpp>
#include<kernel/geometry/conformal_mesh.hpp>
#include<kernel/geometry/cell_sub_set.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;

template<typename Tag_,
         typename IndexType_,
         typename Algo_,
         template<typename, typename> class OT_, typename IT_>
class HaloTest:
  public TaggedTest<Tag_, IndexType_, Algo_>
{
  public:
    HaloTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_, Algo_>("HaloTest<" + tag + ">")
    {
    }

    void run() const
    {

      //##################################################################
      //     0  1
      //   0--1--2     *--*--*
      // 2 | 3|  |4    | 0| 1|
      //   3--4--5     *--*--*
      //    5  6

      Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_> > m3(0);

      //configure attribute
      Foundation::Attribute<double, std::vector> attr;

      //add vertices
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      attr.push_back(double(0));
      attr.push_back(double(0.5));
      attr.push_back(double(1));
      attr.push_back(double(0));
      attr.push_back(double(0.5));
      attr.push_back(double(1));

      //add edges
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);

      //add faces
      m3.add_polytope(Foundation::pl_face);
      m3.add_polytope(Foundation::pl_face);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 0, 0); //v->e is set automagically
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 0, 1);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 1, 1);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 1, 2);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 2, 0);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 2, 3);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 3, 1);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 3, 4);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 4, 2);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 4, 5);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 5, 3);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 5, 4);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 6, 4);
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 6, 5);

      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 0); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 1); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 3); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 4); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 1); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 2); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 4); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 5); //v->f is set automagically

      //clone mesh
      Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_> > m4(1, m3);

      //init simple halo
      Foundation::Halo<0, PLEdge, Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_> > > h(m3, 1);

      //add connections
      //
      // *--*--*
      // |0 | 1| m3
      // *--*--*
      //  5   6
      //  |   |
      //  0   1
      // *--*--*
      // |0 | 1| m4
      // *--*--*

      h.push_back(5u);
      h.push_back(6u);

      TEST_CHECK_EQUAL(h.size(), 2u);
      TEST_CHECK_EQUAL(h.get_element(0u), 5u);
      TEST_CHECK_EQUAL(h.get_element(1u), 6u);

    }
};
HaloTest<Mem::Main, Index, Algo::Generic, std::vector, std::vector<Index> > halo_test_cpu_v_v("std::vector, std::vector");
HaloTest<Mem::Main, Index, Algo::Generic, std::deque, std::vector<Index> > halo_test_cpu_d_v("std::deque, std::vector");
HaloTest<Mem::Main, Index, Algo::Generic, std::vector, std::deque<Index> > halo_test_cpu_v_d("std::vector, std::deque");
HaloTest<Mem::Main, Index, Algo::Generic, std::deque, std::deque<Index> > halo_test_cpu_d_d("std::deque, std::deque");

template<typename Tag_, typename IndexType_, typename Algo_, template<typename, typename> class OT_, typename IT_>
class HaloTestGeometryInterface:
  public TaggedTest<Tag_, IndexType_, Algo_>
{
  public:
    HaloTestGeometryInterface(const std::string & tag) :
      TaggedTest<Tag_, IndexType_, Algo_>("HaloTestGeometryInterface<" + tag + ">")
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


      //creating conformal geometry mesh
      IndexType_ index_set[] = {4, 4, 1};
      typedef Geometry::ConformalMesh<Shape::Hypercube<2> > confmeshtype_;
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
      vertex_coord_tuples[0][0] = ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().at(0); //xcoord of first node
      vertex_coord_tuples[0][1] = ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().at(0); //ycoord of first node
      vertex_coord_tuples[1][0] = ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().at(1);
      vertex_coord_tuples[1][1] = ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().at(1);
      vertex_coord_tuples[2][0] = ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().at(2);
      vertex_coord_tuples[2][1] = ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().at(2);
      vertex_coord_tuples[3][0] = ((Foundation::Attribute<double, OT_>*)(attrs.at(0).get()))->get_data().at(3);
      vertex_coord_tuples[3][1] = ((Foundation::Attribute<double, OT_>*)(attrs.at(1).get()))->get_data().at(3);

      //create halo with one edge-edge pair
      Foundation::Halo<0, PLEdge, Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_> > > h(m, 1);
      h.push_back(3);

      //create CellSubSet
      Index polytopes_in_subset[3] = {2, 1, 0}; //no overlap and one edge means two vertices but no faces
      Geometry::CellSubSet<Shape::Hypercube<2> > cell_sub_set(polytopes_in_subset);

      typename Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_> >::storage_type_ adjacent_vertices(h.get_mesh()->get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, h.get_element(0)));

      cell_sub_set.template get_target_set<0>()[0] = adjacent_vertices.at(0); //first vertex
      cell_sub_set.template get_target_set<0>()[1] = adjacent_vertices.at(1); //second vertex
      cell_sub_set.template get_target_set<1>()[0] = h.get_element(0); //only edge

      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[0], 1ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<0>()[1], 3ul);
      TEST_CHECK_EQUAL(cell_sub_set.template get_target_set<1>()[0], 3ul);
    }
};
HaloTestGeometryInterface<Mem::Main, Index, Algo::Generic, std::vector, std::vector<Index> > halo_test_fginter_cpu_v_v("std::vector, std::vector");
HaloTestGeometryInterface<Mem::Main, Index, Algo::Generic, std::vector, std::deque<Index> > halo_test_fginter_cpu_v_d("std::vector, std::deque");
