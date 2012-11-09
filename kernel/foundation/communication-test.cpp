#include <kernel/base_header.hpp>
#ifdef SERIAL
#include <test_system/test_system.hpp>

#include<kernel/foundation/mesh.hpp>
#include<kernel/foundation/halo.hpp>
#include<kernel/foundation/halo_data.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/foundation/attribute.hpp>
#include<kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class CommunicationTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    CommunicationTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("CommunicationTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      //##################################################################
      //     0  1
      //   0--1--2     *--*--*
      // 2 | 3|  |4    | 0| 1|
      //   3--4--5     *--*--*
      //    5  6

      //build container for all function data associated with m3
      typename Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> >::attr_base_type_ all_attributes_m3;

      //create mesh
      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m3(0, &all_attributes_m3);

      //configure attribute
      all_attributes_m3.push_back(std::shared_ptr<Foundation::AttributeBase<> >(new Foundation::Attribute<double>()));
      Foundation::MeshAttributeRegistration::execute(m3, Foundation::pl_vertex);

      //add vertices
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);

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
      ((Foundation::Attribute<double>*)(all_attributes_m3.at(0).get()))->get_data().push_back(double(42.));
      ((Foundation::Attribute<double>*)(all_attributes_m3.at(0).get()))->get_data().push_back(double(47.));

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
      typename Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> >::attr_base_type_ all_attributes_m4;
      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m4(1, m3, &all_attributes_m4);
      //configure attribute
      all_attributes_m4.push_back(std::shared_ptr<Foundation::AttributeBase<> >(new Foundation::Attribute<double>()));
      Foundation::MeshAttributeRegistration::execute(m4, Foundation::pl_vertex);

      //alter m4's attribute values
      ((Foundation::Attribute<double>*)(all_attributes_m4.at(0).get()))->push_back(3333.);
      ((Foundation::Attribute<double>*)(all_attributes_m4.at(0).get()))->push_back(4444.);

      //init simple halo
      Foundation::Halo<0, Foundation::pl_face, Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > > h(m3, 1);

      //add connections
      //
      // *--*--*
      // |0 | 1| m3
      // *--*--*
      //  |   |
      // *--*--*
      // |0 | 1| m4
      // *--*--*

      h.add_element_pair(0u, 0u);
      h.add_element_pair(1u, 1u);

      TEST_CHECK_EQUAL(h.size(), 2u);
      TEST_CHECK_EQUAL(h.get_element(0u), 0u);
      TEST_CHECK_EQUAL(h.get_element(1u), 1u);
      TEST_CHECK_EQUAL(h.get_element_counterpart(0u), 0u);
      TEST_CHECK_EQUAL(h.get_element_counterpart(1u), 1u);

      Foundation::InterfacedComm<Tag_, Foundation::com_exchange>::execute(h, *(Foundation::Attribute<double>*)(all_attributes_m4.at(0).get()));

      TEST_CHECK_EQUAL(((Foundation::Attribute<double>*)(all_attributes_m3.at(0).get()))->get_data().at(0), 42.);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double>*)(all_attributes_m3.at(0).get()))->get_data().at(1), 47.);
    }
};
CommunicationTest<Archs::CPU, unsigned long, std::vector, std::vector<unsigned long> > halo_test_cpu_v_v("std::vector, std::vector");
CommunicationTest<Archs::CPU, unsigned long, std::deque, std::vector<unsigned long> > halo_test_cpu_d_v("std::deque, std::vector");
CommunicationTest<Archs::CPU, unsigned long, std::vector, std::deque<unsigned long> > halo_test_cpu_v_d("std::vector, std::deque");
CommunicationTest<Archs::CPU, unsigned long, std::deque, std::deque<unsigned long> > halo_test_cpu_d_d("std::deque, std::deque");
#endif // SERIAL
