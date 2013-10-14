#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/partitioning.hpp>
#include<kernel/foundation/load_balancing.hpp>
#include<kernel/archs.hpp>

#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;

template<typename Tag_, typename IndexType_, typename Algo_, template<typename, typename> class OT_, typename IT_>
class LoadBalancingTest1D:
  public TaggedTest<Tag_, IndexType_, Algo_>
{
  public:
    LoadBalancingTest1D(const std::string & tag) :
      TaggedTest<Tag_, Index, Algo_>("LoadBalancingTest1D<" + tag + ">")
    {
    }

    virtual void run() const
    {
      /* (0)  (1)
       *  *----*
       */

      //create attributes for vertex coords
      OT_<Attribute<double, OT_>, std::allocator<Attribute<double, OT_> > > attrs;
      attrs.push_back(Attribute<double, OT_>()); //vertex x-coords

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(0).get_data().push_back(double(1));

      /*
       *  *--0-*
       *  0    1
       */

      //creating foundation mesh
      Foundation::Mesh<Dim1D, Foundation::Topology<IndexType_, OT_, IT_>, OT_> m(0);

      m.add_polytope(Foundation::pl_vertex);
      m.add_polytope(Foundation::pl_vertex);

      m.add_polytope(Foundation::pl_edge);

      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 0, 0);
      m.add_adjacency(Foundation::pl_vertex, Foundation::pl_edge, 1, 0);

      OT_<std::shared_ptr<HaloBase<Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, OT_> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, OT_> > > > halos;
      OT_<Halo<0, PLVertex, Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>, std::allocator<Halo<0, PLVertex, Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, OT_> > > boundaries;

      boundaries.push_back(Halo<0, PLVertex, Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>(m));
      boundaries.push_back(Halo<0, PLVertex, Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>(m));
      boundaries.at(0).push_back(0);
      boundaries.at(1).push_back(1);

      Index num_procs(5);
      Index rank(0);
      //Index level(4);

      PData<Dim1D, Topology<IndexType_, OT_, IT_>, OT_, Mesh, double> p0(Partitioning<Tag_,
                                                                                      Algo_,
                                                                                      Dim1D,
                                                                                      0,
                                                                                      pl_vertex>::execute(m,
                                                                                                          boundaries,
                                                                                                          num_procs, rank,
                                                                                                          attrs
                                                                                                          ));

      LoadBalancing<LBPUniformCompScaledComm<double> >::execute(p0);

      TEST_CHECK_EQUAL_WITHIN_EPS(p0.comm_cost, double(11), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(p0.comp_cost, double(2), std::numeric_limits<double>::epsilon());

    }
};
LoadBalancingTest1D<Mem::Main, Index, Algo::Generic, std::vector, std::vector<Index> > ref_test_cpu_v_v("std::vector, std::vector");
LoadBalancingTest1D<Mem::Main, Index, Algo::Generic, std::vector, std::deque<Index> > ref_test_cpu_v_d("std::vector, std::deque");
LoadBalancingTest1D<Mem::Main, Index, Algo::Generic, std::deque, std::vector<Index> > ref_test_cpu_d_v("std::deque, std::vector");
LoadBalancingTest1D<Mem::Main, Index, Algo::Generic, std::deque, std::deque<Index> > ref_test_cpu_d_d("std::deque, std::deque");

template<typename Tag_, typename IndexType_, typename Algo_, template<typename, typename> class OT_, typename IT_>
class LoadBalancingTest2D:
  public TaggedTest<Tag_, IndexType_, Algo_>
{
  public:
    LoadBalancingTest2D(const std::string & tag) :
      TaggedTest<Tag_, Index, Algo_>("LoadBalancingTest2D<" + tag + ">")
    {
    }

    virtual void run() const
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

      Foundation::Mesh<Dim2D, Foundation::Topology<IndexType_, OT_, IT_>, OT_> m_fine(m);

      OT_<Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>, std::allocator<Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, OT_> > > boundaries;
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>(m));
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>(m));
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>(m));
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, OT_>(m));
      boundaries.at(0).push_back(0);
      boundaries.at(1).push_back(1);
      boundaries.at(2).push_back(2);
      boundaries.at(3).push_back(3);

      Index num_procs(3);
      Index rank(0);
      //Index level(4);

      PData<Dim2D, Topology<IndexType_, OT_, IT_>, OT_, Mesh, double> p0(Partitioning<Tag_,
                                                                                      Algo_,
                                                                                      Dim2D,
                                                                                      0,
                                                                                      pl_vertex>::execute(m,
                                                                                                          boundaries,
                                                                                                          num_procs, rank,
                                                                                                          attrs
                                                                                                          ));

      LoadBalancing<LBPUniformCompScaledComm<double> >::execute(p0);

      TEST_CHECK_EQUAL_WITHIN_EPS(p0.comm_cost, double(48), std::numeric_limits<double>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(p0.comp_cost, double(2), std::numeric_limits<double>::epsilon());
    }
};
LoadBalancingTest2D<Mem::Main, Index, Algo::Generic, std::vector, std::vector<Index> > ref_test_cpu_v_v_2d("std::vector, std::vector");
LoadBalancingTest2D<Mem::Main, Index, Algo::Generic, std::vector, std::deque<Index> > ref_test_cpu_v_d_2d("std::vector, std::deque");
LoadBalancingTest2D<Mem::Main, Index, Algo::Generic, std::deque, std::vector<Index> > ref_test_cpu_d_v_2d("std::deque, std::vector");
LoadBalancingTest2D<Mem::Main, Index, Algo::Generic, std::deque, std::deque<Index> > ref_test_cpu_d_d_2d("std::deque, std::deque");
