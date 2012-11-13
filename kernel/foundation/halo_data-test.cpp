#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/foundation/halo_data.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class HaloDataTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    HaloDataTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("HaloDataTest<" + tag + ">")
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

      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m3(0);

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
      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m4(1, m3);

      //init simple halo
      Foundation::Halo<0, Foundation::pl_edge, Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > > h(m3, 1);

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

      h.add_element_pair(5u, 0u);
      h.add_element_pair(6u, 1u);

      TEST_CHECK_EQUAL(h.size(), 2u);
      TEST_CHECK_EQUAL(h.get_element(0u), 5u);
      TEST_CHECK_EQUAL(h.get_element(1u), 6u);
      TEST_CHECK_EQUAL(h.get_element_counterpart(0u), 0u);
      TEST_CHECK_EQUAL(h.get_element_counterpart(1u), 1u);

      Foundation::HaloData<Foundation::Halo<0,
                                            Foundation::pl_edge,
                                            Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > >,
                           LAFEM::DenseVector,
                           Tag_> hd(h);

    }
};
HaloDataTest<Mem::Main, unsigned long, std::vector, std::vector<unsigned long> > halo_test_cpu_v_v("std::vector, std::vector");
HaloDataTest<Mem::Main, unsigned long, std::deque, std::vector<unsigned long> > halo_test_cpu_d_v("std::deque, std::vector");
HaloDataTest<Mem::Main, unsigned long, std::vector, std::deque<unsigned long> > halo_test_cpu_v_d("std::vector, std::deque");
HaloDataTest<Mem::Main, unsigned long, std::deque, std::deque<unsigned long> > halo_test_cpu_d_d("std::deque, std::deque");
