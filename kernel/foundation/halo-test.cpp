#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<typename Tag_,
         typename IndexType_,
         template<typename, typename> class OT_, typename IT_>
class HaloTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    HaloTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("HaloTest<" + tag + ">")
    {
    }

    void run() const
    {

      std::cout << "FIN" << std::endl;
      //##################################################################
      //     0  1
      //   0--1--2     *--*--*
      // 2 | 3|  |4    | 0| 1|
      //   3--4--5     *--*--*
      //    5  6

      std::cout << "FIN" << std::endl;
      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m3(0);
      std::cout << "FIN" << std::endl;

      //configure attribute
      std::cout << "FIN" << std::endl;
      Foundation::Attribute<double, std::vector> attr;
      Foundation::MeshAttributeRegistration::execute(m3, Foundation::pl_vertex);

      //add vertices
      std::cout << "FIN" << std::endl;
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
      std::cout << "FIN" << std::endl;
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);
      m3.add_polytope(Foundation::pl_edge);

      //add faces
      std::cout << "FIN" << std::endl;
      m3.add_polytope(Foundation::pl_face);
      m3.add_polytope(Foundation::pl_face);

      std::cout << "FIN" << std::endl;
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

      std::cout << "FIN" << std::endl;
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 0); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 1); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 3); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 0, 4); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 1); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 2); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 4); //v->f is set automagically
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_vertex, 1, 5); //v->f is set automagically

      //clone mesh
      std::cout << "FIN" << std::endl;
      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m4(1, m3);
      std::cout << "FIN" << std::endl;

      //init simple halo
      std::cout << "FIN" << std::endl;
      Foundation::Halo<0, Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > > h(m3, 1);
      std::cout << "FIN" << std::endl;

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

      h.add_halo_element_pair(5u, 0u);
      h.add_halo_element_pair(6u, 1u);

      std::cout << "FIN" << std::endl;
      TEST_CHECK_EQUAL(h.size(), 2u);
      std::cout << "FIN" << std::endl;
      TEST_CHECK_EQUAL(h.get_element(0u), 5u);
      std::cout << "FIN" << std::endl;
      TEST_CHECK_EQUAL(h.get_element(1u), 6u);
      std::cout << "FIN" << std::endl;
      TEST_CHECK_EQUAL(h.get_element_counterpart(0u), 0u);
      std::cout << "FIN" << std::endl;
      TEST_CHECK_EQUAL(h.get_element_counterpart(1u), 1u);

      std::cout << "FIN" << std::endl;
    }
};
HaloTest<Archs::None, Index, std::vector, std::vector<Index> > halo_test_cpu_v_v("std::vector, std::vector");
/*HaloTest<Archs::None, Index, std::deque, std::vector<Index> > halo_test_cpu_d_v("std::deque, std::vector");
HaloTest<Archs::None, Index, std::vector, std::deque<Index> > halo_test_cpu_v_d("std::vector, std::deque");
HaloTest<Archs::None, Index, std::deque, std::deque<Index> > halo_test_cpu_d_d("std::deque, std::deque");*/
