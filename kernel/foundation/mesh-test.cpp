#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/mesh.hpp>
#include<kernel/foundation/topology.hpp>
#include<kernel/foundation/dense_data_wrapper.hpp>
#include<kernel/lafem/dense_vector.hpp>
#include<kernel/archs.hpp>
#include<kernel/geometry/conformal_mesh.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MeshTestAttr:
public TaggedTest<Tag_, IndexType_>
{
public:
  MeshTestAttr(const std::string & tag) :
    TaggedTest<Tag_, IndexType_>("MeshTestAttr<" + tag + ">")
  {
  }

  void run() const
  {
    //basic tests
    Foundation::Mesh<> m(0);

    Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m2(1);

    TEST_CHECK_EQUAL(m.get_num_levels(), 3ul);
    TEST_CHECK_EQUAL(m2.get_num_levels(), 3ul);
    //##################################################################

    Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m3(2);

    //configure attribute
    Foundation::Attribute<double, std::vector> attr;
    Foundation::MeshAttributeRegistration::execute(m3, Foundation::pl_vertex);
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

      TEST_CHECK_EQUAL(m3.get_topologies().at(Foundation::ipi_vertex_edge).get_history().get_functors().size(), 6);
      TEST_CHECK_EQUAL(m3.get_topologies().at(Foundation::ipi_vertex_edge).size(), 6u);
      TEST_CHECK_EQUAL(m3.get_topologies().at(Foundation::ipi_vertex_face).size(), 6u);
      TEST_CHECK_EQUAL(m3.get_topologies().at(Foundation::ipi_edge_vertex).size(), 7u);
      TEST_CHECK_EQUAL(m3.get_topologies().at(Foundation::ipi_face_vertex).size(), 2u);


      /*     0  1
           0--1--2     *--*--*
         2 | 3|  |4    | 0| 1|
           3--4--5     *--*--*
            5  6
      */
      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 0, 0); //v->e is set automagically
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 0));
      TEST_CHECK_EQUAL(test_0.size(), 1ul);
      TEST_CHECK_EQUAL(test_0.at(0), 0ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 0));
      TEST_CHECK_EQUAL(test_0_.size(), 1ul);
      TEST_CHECK_EQUAL(test_0_.at(0), 0ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 0, 1);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0a(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 0));
      TEST_CHECK_EQUAL(test_0a.size(), 2ul);
      TEST_CHECK_EQUAL(test_0a.at(0), 0ul);
      TEST_CHECK_EQUAL(test_0a.at(1), 1ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0a_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 1));
      TEST_CHECK_EQUAL(test_0a_.size(), 1ul);
      TEST_CHECK_EQUAL(test_0a_.at(0), 0ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 1, 1);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0b(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 1));
      TEST_CHECK_EQUAL(test_0b.size(), 1ul);
      TEST_CHECK_EQUAL(test_0b.at(0), 1ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0b_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 1));
      TEST_CHECK_EQUAL(test_0b_.size(), 2ul);
      TEST_CHECK_EQUAL(test_0b_.at(0), 0ul);
      TEST_CHECK_EQUAL(test_0b_.at(1), 1ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 1, 2);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0c(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 1));
      TEST_CHECK_EQUAL(test_0c.size(), 2ul);
      TEST_CHECK_EQUAL(test_0c.at(0), 1ul);
      TEST_CHECK_EQUAL(test_0c.at(1), 2ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0c_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 2));
      TEST_CHECK_EQUAL(test_0c_.size(), 1ul);
      TEST_CHECK_EQUAL(test_0c_.at(0), 1ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 2, 0);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0d(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 2));
      TEST_CHECK_EQUAL(test_0d.size(), 1ul);
      TEST_CHECK_EQUAL(test_0d.at(0), 0ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0d_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 0));
      TEST_CHECK_EQUAL(test_0d_.size(), 2ul);
      TEST_CHECK_EQUAL(test_0d_.at(1), 2ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 2, 3);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0e(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 2));
      TEST_CHECK_EQUAL(test_0e.size(), 2ul);
      TEST_CHECK_EQUAL(test_0e.at(0), 0ul);
      TEST_CHECK_EQUAL(test_0e.at(1), 3ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0e_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 3));
      TEST_CHECK_EQUAL(test_0e_.size(), 1ul);
      TEST_CHECK_EQUAL(test_0e_.at(0), 2ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 3, 1);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0f(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 3));
      TEST_CHECK_EQUAL(test_0f.size(), 1ul);
      TEST_CHECK_EQUAL(test_0f.at(0), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0f_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 1));
      TEST_CHECK_EQUAL(test_0f_.size(), 3ul);
      TEST_CHECK_EQUAL(test_0f_.at(2), 3ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 3, 4);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0g(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 3));
      TEST_CHECK_EQUAL(test_0g.size(), 2ul);
      TEST_CHECK_EQUAL(test_0g.at(0), 1ul);
      TEST_CHECK_EQUAL(test_0g.at(1), 4ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0g_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 4));
      TEST_CHECK_EQUAL(test_0g_.size(), 1ul);
      TEST_CHECK_EQUAL(test_0g_.at(0), 3ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 4, 2);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0h(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 4));
      TEST_CHECK_EQUAL(test_0h.size(), 1ul);
      TEST_CHECK_EQUAL(test_0h.at(0), 2ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0h_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 2));
      TEST_CHECK_EQUAL(test_0h_.size(), 2ul);
      TEST_CHECK_EQUAL(test_0h_.at(1), 4ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 4, 5);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0i(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 4));
      TEST_CHECK_EQUAL(test_0i.size(), 2ul);
      TEST_CHECK_EQUAL(test_0i.at(0), 2ul);
      TEST_CHECK_EQUAL(test_0i.at(1), 5ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0i_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 5));
      TEST_CHECK_EQUAL(test_0i_.size(), 1ul);
      TEST_CHECK_EQUAL(test_0i_.at(0), 4ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 5, 3);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0j(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 5));
      TEST_CHECK_EQUAL(test_0j.size(), 1ul);
      TEST_CHECK_EQUAL(test_0j.at(0), 3ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0j_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 3));
      TEST_CHECK_EQUAL(test_0j_.size(), 2ul);
      TEST_CHECK_EQUAL(test_0j_.at(1), 5ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 5, 4);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0k(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 5));
      TEST_CHECK_EQUAL(test_0k.size(), 2ul);
      TEST_CHECK_EQUAL(test_0k.at(0), 3ul);
      TEST_CHECK_EQUAL(test_0k.at(1), 4ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0k_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 4));
      TEST_CHECK_EQUAL(test_0k_.size(), 2ul);
      TEST_CHECK_EQUAL(test_0k_.at(1), 5ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 6, 4);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0l(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 6));
      TEST_CHECK_EQUAL(test_0l.size(), 1ul);
      TEST_CHECK_EQUAL(test_0l.at(0), 4ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0l_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 4));
      TEST_CHECK_EQUAL(test_0l_.size(), 3ul);
      TEST_CHECK_EQUAL(test_0l_.at(2), 6ul);

      m3.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 6, 5);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0m(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 6));
      TEST_CHECK_EQUAL(test_0m.size(), 2ul);
      TEST_CHECK_EQUAL(test_0m.at(0), 4ul);
      TEST_CHECK_EQUAL(test_0m.at(1), 5ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0m_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 5));
      TEST_CHECK_EQUAL(test_0m_.size(), 2ul);
      TEST_CHECK_EQUAL(test_0m_.at(1), 6ul);

//---------------------------------------------------------------------------------------------------
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 0); //v->f is set automagically
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0n(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, 0));
      TEST_CHECK_EQUAL(test_0n.size(), 1ul);
      TEST_CHECK_EQUAL(test_0n.at(0), 0ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0n_(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_face, 0));
      TEST_CHECK_EQUAL(test_0n_.size(), 1ul);
      TEST_CHECK_EQUAL(test_0n_.at(0), 0ul);

      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 1); //v->f is set automagically
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0o(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, 0));
      TEST_CHECK_EQUAL(test_0o.size(), 2ul);
      TEST_CHECK_EQUAL(test_0o.at(1), 1ul);

      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 3); //v->f is set automagically
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0p(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, 0));
      TEST_CHECK_EQUAL(test_0p.size(), 3ul);
      TEST_CHECK_EQUAL(test_0p.at(2), 3ul);

      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 4); //v->f is set automagically
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0q(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, 0));
      TEST_CHECK_EQUAL(test_0q.size(), 4ul);
      TEST_CHECK_EQUAL(test_0q.at(0), 0ul);
      TEST_CHECK_EQUAL(test_0q.at(1), 1ul);
      TEST_CHECK_EQUAL(test_0q.at(2), 3ul);
      TEST_CHECK_EQUAL(test_0q.at(3), 4ul);

      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 1); //v->f is set automagically
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0r(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, 1));
      TEST_CHECK_EQUAL(test_0r.size(), 1ul);
      TEST_CHECK_EQUAL(test_0r.at(0), 1ul);

      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 2); //v->f is set automagically
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0s(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, 1));
      TEST_CHECK_EQUAL(test_0s.size(), 2ul);
      TEST_CHECK_EQUAL(test_0s.at(1), 2ul);

      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 4); //v->f is set automagically
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0t(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, 1));
      TEST_CHECK_EQUAL(test_0t.size(), 3ul);
      TEST_CHECK_EQUAL(test_0t.at(2), 4ul);

      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 5); //v->f is set automagically
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0u(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, 1));
      TEST_CHECK_EQUAL(test_0u.size(), 4ul);
      TEST_CHECK_EQUAL(test_0u.at(0), 1ul);
      TEST_CHECK_EQUAL(test_0u.at(1), 2ul);
      TEST_CHECK_EQUAL(test_0u.at(2), 4ul);
      TEST_CHECK_EQUAL(test_0u.at(3), 5ul);


      //testing face-edge access
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_1(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_edge, 0));
      TEST_CHECK_EQUAL(test_1.size(), 6ul);
      TEST_CHECK_EQUAL(test_1.at(0), 0ul);
      TEST_CHECK_EQUAL(test_1.at(1), 2ul);
      TEST_CHECK_EQUAL(test_1.at(2), 1ul);
      TEST_CHECK_EQUAL(test_1.at(3), 3ul);
      TEST_CHECK_EQUAL(test_1.at(4), 5ul);
      TEST_CHECK_EQUAL(test_1.at(5), 6ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_2(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_edge, 1));
      TEST_CHECK_EQUAL(test_2.size(), 6ul);
      TEST_CHECK_EQUAL(test_2.at(0), 0ul);
      TEST_CHECK_EQUAL(test_2.at(1), 1ul);
      TEST_CHECK_EQUAL(test_2.at(2), 3ul);
      TEST_CHECK_EQUAL(test_2.at(3), 4ul);
      TEST_CHECK_EQUAL(test_2.at(4), 5ul);
      TEST_CHECK_EQUAL(test_2.at(5), 6ul);

      //testing face-face access
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_3(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_face, 0));
      TEST_CHECK_EQUAL(test_3.size(), 2ul);
      TEST_CHECK_EQUAL(test_3.at(0), 0ul);
      TEST_CHECK_EQUAL(test_3.at(1), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_4(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_face, 1));
      TEST_CHECK_EQUAL(test_4.size(), 2ul);
      TEST_CHECK_EQUAL(test_4.at(0), 0ul);
      TEST_CHECK_EQUAL(test_4.at(1), 1ul);

      //testing face-vertex access
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_5(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, 0));
      TEST_CHECK_EQUAL(test_5.size(), 4ul);
      TEST_CHECK_EQUAL(test_5.at(0), 0ul);
      TEST_CHECK_EQUAL(test_5.at(1), 1ul);
      TEST_CHECK_EQUAL(test_5.at(2), 3ul);
      TEST_CHECK_EQUAL(test_5.at(3), 4ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_6(m3.get_adjacent_polytopes(Foundation::pl_face, Foundation::pl_vertex, 1));
      TEST_CHECK_EQUAL(test_6.size(), 4ul);
      TEST_CHECK_EQUAL(test_6.at(0), 1ul);
      TEST_CHECK_EQUAL(test_6.at(1), 2ul);
      TEST_CHECK_EQUAL(test_6.at(2), 4ul);
      TEST_CHECK_EQUAL(test_6.at(3), 5ul);

      //testing edge-vertex
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_7(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 0));
      TEST_CHECK_EQUAL(test_7.size(), 2ul);
      TEST_CHECK_EQUAL(test_7.at(0), 0ul);
      TEST_CHECK_EQUAL(test_7.at(1), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_8(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 1));
      TEST_CHECK_EQUAL(test_8.size(), 2ul);
      TEST_CHECK_EQUAL(test_8.at(0), 1ul);
      TEST_CHECK_EQUAL(test_8.at(1), 2ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_9(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 2));
      TEST_CHECK_EQUAL(test_9.size(), 2ul);
      TEST_CHECK_EQUAL(test_9.at(0), 0ul);
      TEST_CHECK_EQUAL(test_9.at(1), 3ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_10(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 3));
      TEST_CHECK_EQUAL(test_10.size(), 2ul);
      TEST_CHECK_EQUAL(test_10.at(0), 1ul);
      TEST_CHECK_EQUAL(test_10.at(1), 4ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_11(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 4));
      TEST_CHECK_EQUAL(test_11.size(), 2ul);
      TEST_CHECK_EQUAL(test_11.at(0), 2ul);
      TEST_CHECK_EQUAL(test_11.at(1), 5ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_12(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 5));
      TEST_CHECK_EQUAL(test_12.size(), 2ul);
      TEST_CHECK_EQUAL(test_12.at(0), 3ul);
      TEST_CHECK_EQUAL(test_12.at(1), 4ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_13(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 6));
      TEST_CHECK_EQUAL(test_13.size(), 2ul);
      TEST_CHECK_EQUAL(test_13.at(0), 4ul);
      TEST_CHECK_EQUAL(test_13.at(1), 5ul);

      //testing edge-edge
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_14(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_edge, 0));
      TEST_CHECK_EQUAL(test_14.size(), 4ul);
      TEST_CHECK_EQUAL(test_14.at(0), 0ul);
      TEST_CHECK_EQUAL(test_14.at(1), 2ul);
      TEST_CHECK_EQUAL(test_14.at(2), 1ul);
      TEST_CHECK_EQUAL(test_14.at(3), 3ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_15(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_edge, 1));
      TEST_CHECK_EQUAL(test_15.size(), 4ul);
      TEST_CHECK_EQUAL(test_15.at(0), 0ul);
      TEST_CHECK_EQUAL(test_15.at(1), 1ul);
      TEST_CHECK_EQUAL(test_15.at(2), 3ul);
      TEST_CHECK_EQUAL(test_15.at(3), 4ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_16(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_edge, 2));
      TEST_CHECK_EQUAL(test_16.size(), 3ul);
      TEST_CHECK_EQUAL(test_16.at(0), 0ul);
      TEST_CHECK_EQUAL(test_16.at(1), 2ul);
      TEST_CHECK_EQUAL(test_16.at(2), 5ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_17(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_edge, 3));
      TEST_CHECK_EQUAL(test_17.size(), 5ul);
      TEST_CHECK_EQUAL(test_17.at(0), 0ul);
      TEST_CHECK_EQUAL(test_17.at(1), 1ul);
      TEST_CHECK_EQUAL(test_17.at(2), 3ul);
      TEST_CHECK_EQUAL(test_17.at(3), 5ul);
      TEST_CHECK_EQUAL(test_17.at(4), 6ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_18(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_edge, 4));
      TEST_CHECK_EQUAL(test_18.size(), 3ul);
      TEST_CHECK_EQUAL(test_18.at(0), 1ul);
      TEST_CHECK_EQUAL(test_18.at(1), 4ul);
      TEST_CHECK_EQUAL(test_18.at(2), 6ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_19(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_edge, 5));
      TEST_CHECK_EQUAL(test_19.size(), 4ul);
      TEST_CHECK_EQUAL(test_19.at(0), 2ul);
      TEST_CHECK_EQUAL(test_19.at(1), 5ul);
      TEST_CHECK_EQUAL(test_19.at(2), 3ul);
      TEST_CHECK_EQUAL(test_19.at(3), 6ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_20(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_edge, 6));
      TEST_CHECK_EQUAL(test_20.size(), 4ul);
      TEST_CHECK_EQUAL(test_20.at(0), 3ul);
      TEST_CHECK_EQUAL(test_20.at(1), 5ul);
      TEST_CHECK_EQUAL(test_20.at(2), 6ul);
      TEST_CHECK_EQUAL(test_20.at(3), 4ul);

      //testing edge-face access
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_21(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_face, 0));
      TEST_CHECK_EQUAL(test_21.size(), 2ul);
      TEST_CHECK_EQUAL(test_21.at(0), 0ul);
      TEST_CHECK_EQUAL(test_21.at(1), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_22(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_face, 1));
      TEST_CHECK_EQUAL(test_22.size(), 2ul);
      TEST_CHECK_EQUAL(test_22.at(0), 0ul);
      TEST_CHECK_EQUAL(test_22.at(1), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_23(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_face, 2));
      TEST_CHECK_EQUAL(test_23.size(), 1ul);
      TEST_CHECK_EQUAL(test_23.at(0), 0ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_24(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_face, 3));
      TEST_CHECK_EQUAL(test_24.size(), 2ul);
      TEST_CHECK_EQUAL(test_24.at(0), 0ul);
      TEST_CHECK_EQUAL(test_24.at(1), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_25(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_face, 4));
      TEST_CHECK_EQUAL(test_25.size(), 1ul);
      TEST_CHECK_EQUAL(test_25.at(0), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_26(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_face, 5));
      TEST_CHECK_EQUAL(test_26.size(), 2ul);
      TEST_CHECK_EQUAL(test_26.at(0), 0ul);
      TEST_CHECK_EQUAL(test_26.at(1), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_27(m3.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_face, 6));
      TEST_CHECK_EQUAL(test_27.size(), 2ul);
      TEST_CHECK_EQUAL(test_27.at(0), 0ul);
      TEST_CHECK_EQUAL(test_27.at(1), 1ul);

      //testing vertex-vertex access
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_28(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_vertex, 0));
      TEST_CHECK_EQUAL(test_28.size(), 3ul);
      TEST_CHECK_EQUAL(test_28.at(0), 0ul);
      TEST_CHECK_EQUAL(test_28.at(1), 1ul);
      TEST_CHECK_EQUAL(test_28.at(2), 3ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_29(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_vertex, 1));
      TEST_CHECK_EQUAL(test_29.size(), 4ul);
      TEST_CHECK_EQUAL(test_29.at(0), 0ul);
      TEST_CHECK_EQUAL(test_29.at(1), 1ul);
      TEST_CHECK_EQUAL(test_29.at(2), 2ul);
      TEST_CHECK_EQUAL(test_29.at(3), 4ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_30(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_vertex, 2));
      TEST_CHECK_EQUAL(test_30.size(), 3ul);
      TEST_CHECK_EQUAL(test_30.at(0), 1ul);
      TEST_CHECK_EQUAL(test_30.at(1), 2ul);
      TEST_CHECK_EQUAL(test_30.at(2), 5ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_31(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_vertex, 3));
      TEST_CHECK_EQUAL(test_31.size(), 3ul);
      TEST_CHECK_EQUAL(test_31.at(0), 0ul);
      TEST_CHECK_EQUAL(test_31.at(1), 3ul);
      TEST_CHECK_EQUAL(test_31.at(2), 4ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_32(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_vertex, 4));
      TEST_CHECK_EQUAL(test_32.size(), 4ul);
      TEST_CHECK_EQUAL(test_32.at(0), 1ul);
      TEST_CHECK_EQUAL(test_32.at(1), 4ul);
      TEST_CHECK_EQUAL(test_32.at(2), 3ul);
      TEST_CHECK_EQUAL(test_32.at(3), 5ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_33(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_vertex, 5));
      TEST_CHECK_EQUAL(test_33.size(), 3ul);
      TEST_CHECK_EQUAL(test_33.at(0), 2ul);
      TEST_CHECK_EQUAL(test_33.at(1), 5ul);
      TEST_CHECK_EQUAL(test_33.at(2), 4ul);

      //testing vertex-edge access
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_34(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 0));
      TEST_CHECK_EQUAL(test_34.size(), 2ul);
      TEST_CHECK_EQUAL(test_34.at(0), 0ul);
      TEST_CHECK_EQUAL(test_34.at(1), 2ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_35(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 1));
      TEST_CHECK_EQUAL(test_35.size(), 3ul);
      TEST_CHECK_EQUAL(test_35.at(0), 0ul);
      TEST_CHECK_EQUAL(test_35.at(1), 1ul);
      TEST_CHECK_EQUAL(test_35.at(2), 3ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_36(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 2));
      TEST_CHECK_EQUAL(test_36.size(), 2ul);
      TEST_CHECK_EQUAL(test_36.at(0), 1ul);
      TEST_CHECK_EQUAL(test_36.at(1), 4ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_37(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 3));
      TEST_CHECK_EQUAL(test_37.size(), 2ul);
      TEST_CHECK_EQUAL(test_37.at(0), 2ul);
      TEST_CHECK_EQUAL(test_37.at(1), 5ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_38(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 4));
      TEST_CHECK_EQUAL(test_38.size(), 3ul);
      TEST_CHECK_EQUAL(test_38.at(0), 3ul);
      TEST_CHECK_EQUAL(test_38.at(1), 5ul);
      TEST_CHECK_EQUAL(test_38.at(2), 6ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_39(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 5));
      TEST_CHECK_EQUAL(test_39.size(), 2ul);
      TEST_CHECK_EQUAL(test_39.at(0), 4ul);
      TEST_CHECK_EQUAL(test_39.at(1), 6ul);

      //testing vertex-face access
      try
      {
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_40(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_face, 0));
      TEST_CHECK_EQUAL(test_40.size(), 1ul);
      TEST_CHECK_EQUAL(test_40.at(0), 0ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_41(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_face, 1));
      TEST_CHECK_EQUAL(test_41.size(), 2ul);
      TEST_CHECK_EQUAL(test_41.at(0), 0ul);
      TEST_CHECK_EQUAL(test_41.at(1), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_42(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_face, 2));
      TEST_CHECK_EQUAL(test_42.size(), 1ul);
      TEST_CHECK_EQUAL(test_42.at(0), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_43(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_face, 3));
      TEST_CHECK_EQUAL(test_43.size(), 1ul);
      TEST_CHECK_EQUAL(test_43.at(0), 0ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_44(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_face, 4));
      TEST_CHECK_EQUAL(test_44.size(), 2ul);
      TEST_CHECK_EQUAL(test_44.at(0), 0ul);
      TEST_CHECK_EQUAL(test_44.at(1), 1ul);

      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_45(m3.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_face, 5));
      TEST_CHECK_EQUAL(test_45.size(), 1ul);
      TEST_CHECK_EQUAL(test_45.at(0), 1ul);
      }
      catch(std::exception e)
      {
        std::cout << "testing v f" << std::endl;
      }

      //testing copy-ctor
      try
      {
      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m4(3, m3);
      }
      catch(std::exception e)
      {
        std::cout << "copy" << std::endl;
      }
    }
};
MeshTestAttr<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > mesh_test_cpu_v_v("std::vector, std::vector");
#ifdef DEBUG
MeshTestAttr<Archs::None, unsigned long, std::deque, std::vector<unsigned long> > mesh_test_cpu_d_v("std::deque, std::vector");
MeshTestAttr<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > mesh_test_cpu_v_d("std::vector, std::deque");
MeshTestAttr<Archs::None, unsigned long, std::deque, std::deque<unsigned long> > mesh_test_cpu_d_d("std::deque, std::deque");
#endif
//MeshTestAttr<Archs::CPU, unsigned long, std::vector, Foundation::DenseDataWrapper<100, Archs::CPU, unsigned long, DenseVector> > mesh_test_cpu_v_ddw_DV("std::vector, hornet::DV"); //TODO: ddw needs erase(i) member

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MeshTestHistory:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    MeshTestHistory(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("MeshTestHistory<" + tag + ">")
    {
    }

    void run() const
    {
      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m(0);

      //add vertices
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);
      //add edges
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);
      m.add_polytope(Foundation::pl_edge);

      //add faces
      m.add_polytope(Foundation::pl_face);
      m.add_polytope(Foundation::pl_face);

      /*     0  1
           0--1--2     *--*--*
         2 | 3|  |4    | 0| 1|
           3--4--5     *--*--*
            5  6
      */

      m.add_adjacency(Foundation::pl_edge, Foundation::pl_vertex, 0, 0); //v->e is set automagically
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0(m.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 0));
      TEST_CHECK_EQUAL(test_0.size(), 1ul);
      TEST_CHECK_EQUAL(test_0.at(0), 0ul);
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_0_(m.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 0));
      TEST_CHECK_EQUAL(test_0_.size(), 1ul);
      TEST_CHECK_EQUAL(test_0_.at(0), 0ul);

      //undo the adjacency
      m.undo();
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_1(m.get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, 0));
      typename Foundation::Topology<IndexType_, OT_, IT_>::storage_type_ test_1_(m.get_adjacent_polytopes(Foundation::pl_vertex, Foundation::pl_edge, 0));
      TEST_CHECK_EQUAL(test_1.size(), 0ul);
      TEST_CHECK_EQUAL(test_1_.size(), 0ul);

      //clear all manually
      Foundation::SmartPointer<Foundation::FunctorBase> f1(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f2(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f3(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f4(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f5(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f6(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f7(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f8(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f9(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f10(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f11(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f12(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f13(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f14(m.undo());
      Foundation::SmartPointer<Foundation::FunctorBase> f15(m.undo());

      TEST_CHECK_EQUAL(m.get_history().get_functors().size(), 0);

      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_vertex_edge).get_topology().size(), 0ul);
      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_edge_vertex).get_topology().size(), 0ul);
      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_face_vertex).get_topology().size(), 0ul);
      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_vertex_face).get_topology().size(), 0ul);

      //redo all manually
      m.redo(f1);
      m.redo(f2);
      m.redo(f3);
      m.redo(f4);
      m.redo(f5);
      m.redo(f6);
      m.redo(f7);
      m.redo(f8);
      m.redo(f9);
      m.redo(f10);
      m.redo(f11);
      m.redo(f12);
      m.redo(f13);
      m.redo(f14);
      m.redo(f15);

      TEST_CHECK_EQUAL(m.get_history().get_functors().size(), 15);
      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_vertex_edge).get_topology().size(), 6ul);
      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_vertex_face).get_topology().size(), 6ul);
      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_edge_vertex).get_topology().size(), 7ul);
      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_face_vertex).get_topology().size(), 2ul);

      //clear all again
      m.clear();

      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_vertex_edge).get_topology().size(), 0ul);
      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_vertex_face).get_topology().size(), 0ul);
      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_edge_vertex).get_topology().size(), 0ul);
      TEST_CHECK_EQUAL(m.get_topologies().at(Foundation::ipi_face_vertex).get_topology().size(), 0ul);
    }
};
MeshTestHistory<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > mesh_test_his_cpu_v_v("std::vector, std::vector");

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MeshTestGeometryInterface:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    MeshTestGeometryInterface(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("MeshTestGeometryInterface<" + tag + ">")
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
      std::vector<Foundation::SmartPointer<Foundation::AttributeBase<OT_> > > attrs;
      attrs.push_back(Foundation::SmartPointer<Foundation::AttributeBase<OT_> >(new Foundation::Attribute<double, OT_>)); //vertex x-coords
      attrs.push_back(Foundation::SmartPointer<Foundation::AttributeBase<OT_> >(new Foundation::Attribute<double, OT_>)); //vertex y-coords

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
      typedef Geometry::ConformalMesh<Geometry::ConformalMeshPolicy<Geometry::Shape::Hypercube<2> > > confmeshtype_;
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

      typename confmeshtype_::VertexSetType& vertex_coord_tuples(geo_m.get_vertex_set());
      vertex_coord_tuples[0][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(0); //xcoord of first node
      vertex_coord_tuples[0][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(0); //ycoord of first node
      vertex_coord_tuples[1][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(1);
      vertex_coord_tuples[1][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(1);
      vertex_coord_tuples[2][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(2);
      vertex_coord_tuples[2][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(2);
      vertex_coord_tuples[3][0] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(0).get()))->get_data().at(3);
      vertex_coord_tuples[3][1] = ((Foundation::Attribute<double, OT_>*)(m.get_attributes()->at(1).get()))->get_data().at(3);

      TEST_CHECK_EQUAL(vertex_coord_tuples[0][0], 0);
      TEST_CHECK_EQUAL(vertex_coord_tuples[0][1], 0);
      TEST_CHECK_EQUAL(vertex_coord_tuples[1][0], 1);
      TEST_CHECK_EQUAL(vertex_coord_tuples[1][1], 0);
      TEST_CHECK_EQUAL(vertex_coord_tuples[2][0], 0);
      TEST_CHECK_EQUAL(vertex_coord_tuples[2][1], 1);
      TEST_CHECK_EQUAL(vertex_coord_tuples[3][0], 1);
      TEST_CHECK_EQUAL(vertex_coord_tuples[3][1], 1);
    }
};
MeshTestGeometryInterface<Archs::None, unsigned long, std::vector, std::vector<unsigned long> > mesh_test_fginter_cpu_v_v("std::vector, std::vector");
MeshTestGeometryInterface<Archs::None, unsigned long, std::vector, std::deque<unsigned long> > mesh_test_fginter_cpu_v_d("std::vector, std::deque");
