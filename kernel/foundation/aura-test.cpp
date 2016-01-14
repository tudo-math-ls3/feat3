#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/foundation/base.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/aura.hpp>
#include <kernel/foundation/attribute.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;

template<typename Tag_,
         typename IndexType_,
         template<typename, typename> class OT_, typename IT_>
class AuraTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    AuraTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("AuraTest<" + tag + ">")
    {
    }

    void run() const override
    {
      /*(0,1) (1,1)
       *  *----*
       *  |    |
       *  |    |
       *  *----*
       *(0,0) (1,0)
       */

      //create attributes for vertex coords
      OT_<Attribute<double, OT_>, std::allocator<Attribute<double, OT_> > > attrs;
      attrs.push_back(Attribute<double, OT_>()); //vertex x-coords
      attrs.push_back(Attribute<double, OT_>()); //vertex y-coords

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(1));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(1));

      /*  2    3
       *  *-1--*
       *  2    |
       *  |    3
       *  *--0-*
       *  0    1
       */

      //creating foundation mesh
      Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_>, OT_> m(0);
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

      Halo<0, PLEdge, Mesh<Dim2D> > h0(m);
      decltype(h0) h1(m);
      Halo<0, PLVertex, Mesh<Dim2D> > h2(m);
      decltype(h2) h3(m);

      h0.push_back(0);
      h1.push_back(1);
      h2.push_back(1);
      h2.push_back(3);
      h3.push_back(0);
      h3.push_back(2);

      std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D> > > > halos;
      halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim2D> > >(new Halo<0, PLEdge, Mesh<Dim2D> >(h0)));
      halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim2D> > >(new Halo<0, PLEdge, Mesh<Dim2D> >(h1)));
      halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim2D> > >(new Halo<0, PLVertex, Mesh<Dim2D> >(h2)));
      halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim2D> > >(new Halo<0, PLVertex, Mesh<Dim2D> >(h3)));

      auto h_new(Aura<Tag_, decltype(h0)>::value(halos));

      TEST_CHECK(std::find(h_new.get_elements().begin(), h_new.get_elements().end(), 0) !=  h_new.get_elements().end());
      TEST_CHECK(std::find(h_new.get_elements().begin(), h_new.get_elements().end(), 1) !=  h_new.get_elements().end());
    }
};
AuraTest<Mem::Main, Index, std::vector, std::vector<Index> > aura_test_cpu_v_v("std::vector, std::vector");
