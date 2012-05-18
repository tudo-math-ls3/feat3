#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/mesh_data.hpp>
#include <kernel/foundation/dense_data_wrapper.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;

//Test container
template<typename DT_>
class TestArrayClass
{
  public:
    TestArrayClass(Index size) :
      _size(size),
      _data(new DT_[size])
    {
    }

    DT_ & operator[] (Index i)
    {
      return _data[i];
    }

    Index size()
    {
      return _size;
    }

  private:
    Index _size;
    DT_ * _data;

};

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MeshDataTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    MeshDataTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("MeshDataTest<" + tag + ">")
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

      Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> > m3(0);

      //configure attribute
      unsigned my_attribute_index(Foundation::MeshAttributeRegistration<Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> >, double>::execute(m3, Foundation::pl_vertex));

      //add vertices
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_polytope(Foundation::pl_vertex);
      m3.add_attribute_value(my_attribute_index, double(0));
      m3.add_attribute_value(my_attribute_index, double(0.5));
      m3.add_attribute_value(my_attribute_index, double(1));
      m3.add_attribute_value(my_attribute_index, double(0));
      m3.add_attribute_value(my_attribute_index, double(0.5));
      m3.add_attribute_value(my_attribute_index, double(1));

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

      m3.add_adjacency(Foundation::pl_face, Foundation::pl_edge, 0, 0);
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_edge, 0, 2);
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_edge, 0, 3);
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_edge, 0, 5);
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_edge, 1, 1);
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_edge, 1, 3);
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_edge, 1, 4);
      m3.add_adjacency(Foundation::pl_face, Foundation::pl_edge, 1, 6);

      Foundation::MeshData<Foundation::Mesh<Foundation::rnt_2D, Foundation::Topology<IndexType_, OT_, IT_> >, TestArrayClass > data(m3);
      TestArrayClass<double> vertex_attr = data.get_attribute_of_type_1(0);
      TEST_CHECK_EQUAL(vertex_attr.size(), 6u);
    }
};
MeshDataTest<Nil, unsigned long, std::vector, std::vector<unsigned long> > meshdata_test_cpu_v_v("std::vector, std::vector");
MeshDataTest<Nil, unsigned long, std::deque, std::vector<unsigned long> > meshdata_test_cpu_d_v("std::deque, std::vector");
MeshDataTest<Nil, unsigned long, std::vector, std::deque<unsigned long> > meshdata_test_cpu_v_d("std::vector, std::deque");
MeshDataTest<Nil, unsigned long, std::deque, std::deque<unsigned long> > meshdata_test_cpu_d_d("std::deque, std::deque");
