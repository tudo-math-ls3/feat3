#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/partitioning.hpp>
#include<kernel/archs.hpp>

#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;

template<typename Tag_, typename IndexType_, typename Algo_, template<typename, typename> class OT_, typename IT_>
class PartitioningTest1D:
  public TaggedTest<Tag_, IndexType_, Algo_>
{
  public:
    PartitioningTest1D(const std::string & tag) :
      TaggedTest<Tag_, Index, Algo_>("PartitioningTest1D<" + tag + ">")
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
      //inspect p0 for (1)submesh data, (2)comm_halos, (3)boundaries, (4)attributes
      //(1)submesh data on process 0
      TEST_CHECK_EQUAL(p0.submesh->get_topologies().at(ipi_edge_vertex).size(), 2ul);
      TEST_CHECK_EQUAL(p0.submesh->get_topologies().at(ipi_vertex_edge).size(), 3ul);
      TEST_CHECK_EQUAL(p0.submesh->get_map().size(), 2ul);
      TEST_CHECK_EQUAL(p0.submesh->get_map().at(0), 0ul);
      TEST_CHECK_EQUAL(p0.submesh->get_map().at(1), 4ul);
      auto sm_edge_0_adj_vert(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, 0));
      TEST_CHECK_EQUAL(sm_edge_0_adj_vert.size(), 2ul);
      TEST_CHECK_EQUAL(sm_edge_0_adj_vert.at(0), 0ul);
      TEST_CHECK_EQUAL(sm_edge_0_adj_vert.at(1), 1ul);
      auto sm_edge_1_adj_vert(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, 1));
      TEST_CHECK_EQUAL(sm_edge_1_adj_vert.size(), 2ul);
      TEST_CHECK_EQUAL(sm_edge_1_adj_vert.at(0), 1ul);
      TEST_CHECK_EQUAL(sm_edge_1_adj_vert.at(1), 2ul);

      //(2)comm halos for process 0 (should be one vertex halo<0>)
      TEST_CHECK_EQUAL(p0.comm_halos.size(), 1ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(0)->get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(0)->get_level(), pl_vertex);
      TEST_CHECK_EQUAL(p0.comm_halos.at(0)->size(), 1ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(0)->get_element(0), 2ul);

      //(3)boundaries for process 0 (should be one vertex Halo<0>)
      TEST_CHECK_EQUAL(p0.boundaries.size(), 1ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_level(), pl_vertex);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).size(), 1ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_element(0), 0ul);

      //(4)attributes
      TEST_CHECK_EQUAL(p0.attrs.at(0)->size(), 3ul);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(0).get())->at(0), attrs.at(0).at(0));
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(0).get())->at(1), attrs.at(0).at(5));
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(0).get())->at(2), attrs.at(0).at(3));

    }
};
PartitioningTest1D<Mem::Main, Index, Algo::Generic, std::vector, std::vector<Index> > ref_test_cpu_v_v("std::vector, std::vector");
PartitioningTest1D<Mem::Main, Index, Algo::Generic, std::vector, std::deque<Index> > ref_test_cpu_v_d("std::vector, std::deque");
PartitioningTest1D<Mem::Main, Index, Algo::Generic, std::deque, std::vector<Index> > ref_test_cpu_d_v("std::deque, std::vector");
PartitioningTest1D<Mem::Main, Index, Algo::Generic, std::deque, std::deque<Index> > ref_test_cpu_d_d("std::deque, std::deque");

template<typename Tag_, typename IndexType_, typename Algo_, template<typename, typename> class OT_, typename IT_>
class PartitioningTest2D:
  public TaggedTest<Tag_, IndexType_, Algo_>
{
  public:
    PartitioningTest2D(const std::string & tag) :
      TaggedTest<Tag_, Index, Algo_>("PartitioningTest2D<" + tag + ">")
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

      //inspect p0 for (1)submesh data, (2)comm_halos, (3)boundaries, (4)attributes
      //(1)submesh data on process 0
      TEST_CHECK_EQUAL(p0.submesh->get_topologies().at(ipi_face_vertex).size(), 2ul);
      TEST_CHECK_EQUAL(p0.submesh->get_topologies().at(ipi_edge_vertex).size(), 7ul);
      TEST_CHECK_EQUAL(p0.submesh->get_topologies().at(ipi_vertex_edge).size(), 6ul);

      TEST_CHECK_EQUAL(p0.submesh->get_map().size(), 2ul);
      TEST_CHECK_EQUAL(p0.submesh->get_map().at(0), 0ul);
      TEST_CHECK_EQUAL(p0.submesh->get_map().at(1), 1ul);

      auto sm_face_0_adj_vert(p0.submesh->get_adjacent_polytopes(pl_face, pl_vertex, 0));
      TEST_CHECK_EQUAL(sm_face_0_adj_vert.size(), 4ul);
      TEST_CHECK_EQUAL(sm_face_0_adj_vert.at(0), 0ul);
      TEST_CHECK_EQUAL(sm_face_0_adj_vert.at(1), 1ul);
      TEST_CHECK_EQUAL(sm_face_0_adj_vert.at(2), 2ul);
      TEST_CHECK_EQUAL(sm_face_0_adj_vert.at(3), 5ul);

      auto sm_face_1_adj_vert(p0.submesh->get_adjacent_polytopes(pl_face, pl_vertex, 1));
      TEST_CHECK_EQUAL(sm_face_1_adj_vert.size(), 4ul);
      TEST_CHECK_EQUAL(sm_face_1_adj_vert.at(0), 1ul);
      TEST_CHECK_EQUAL(sm_face_1_adj_vert.at(1), 3ul);
      TEST_CHECK_EQUAL(sm_face_1_adj_vert.at(2), 5ul);
      TEST_CHECK_EQUAL(sm_face_1_adj_vert.at(3), 4ul);

      //intersection in {1,5} implies comm_intersection on vertex level must be {1,5}
      auto vertex_comm_int_for_faces_0_1(p0.submesh->get_comm_intersection(pl_face, pl_vertex, 0, 1));
      TEST_CHECK_EQUAL(vertex_comm_int_for_faces_0_1.size(), 2ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_faces_0_1.at(0), 1ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_faces_0_1.at(1), 5ul);

      auto sm_edge_0_adj_vert(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, 0));
      TEST_CHECK_EQUAL(sm_edge_0_adj_vert.size(), 2ul);
      TEST_CHECK_EQUAL(sm_edge_0_adj_vert.at(0), 0ul);
      TEST_CHECK_EQUAL(sm_edge_0_adj_vert.at(1), 1ul);

      auto sm_edge_1_adj_vert(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, 1));
      TEST_CHECK_EQUAL(sm_edge_1_adj_vert.size(), 2ul);
      TEST_CHECK_EQUAL(sm_edge_1_adj_vert.at(0), 0ul);
      TEST_CHECK_EQUAL(sm_edge_1_adj_vert.at(1), 2ul);

      auto sm_edge_2_adj_vert(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, 2));
      TEST_CHECK_EQUAL(sm_edge_2_adj_vert.size(), 2ul);
      TEST_CHECK_EQUAL(sm_edge_2_adj_vert.at(0), 3ul);
      TEST_CHECK_EQUAL(sm_edge_2_adj_vert.at(1), 4ul);

      auto sm_edge_3_adj_vert(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, 3));
      TEST_CHECK_EQUAL(sm_edge_3_adj_vert.size(), 2ul);
      TEST_CHECK_EQUAL(sm_edge_3_adj_vert.at(0), 1ul);
      TEST_CHECK_EQUAL(sm_edge_3_adj_vert.at(1), 3ul);

      auto sm_edge_4_adj_vert(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, 4));
      TEST_CHECK_EQUAL(sm_edge_4_adj_vert.size(), 2ul);
      TEST_CHECK_EQUAL(sm_edge_4_adj_vert.at(0), 1ul);
      TEST_CHECK_EQUAL(sm_edge_4_adj_vert.at(1), 5ul);

      auto sm_edge_5_adj_vert(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, 5));
      TEST_CHECK_EQUAL(sm_edge_5_adj_vert.size(), 2ul);
      TEST_CHECK_EQUAL(sm_edge_5_adj_vert.at(0), 2ul);
      TEST_CHECK_EQUAL(sm_edge_5_adj_vert.at(1), 5ul);

      auto sm_edge_6_adj_vert(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, 6));
      TEST_CHECK_EQUAL(sm_edge_6_adj_vert.size(), 2ul);
      TEST_CHECK_EQUAL(sm_edge_6_adj_vert.at(0), 4ul);
      TEST_CHECK_EQUAL(sm_edge_6_adj_vert.at(1), 5ul);

      //=> comm_intersection on edge level must be {4}
      auto edge_comm_int_for_faces_0_1(p0.submesh->get_comm_intersection(pl_face, pl_edge, 0, 1));
      TEST_CHECK_EQUAL(edge_comm_int_for_faces_0_1.size(), 1ul);
      TEST_CHECK_EQUAL(edge_comm_int_for_faces_0_1.at(0), 4ul);

      //(2)comm halos for process 0
      TEST_CHECK_EQUAL(p0.comm_halos.size(), 4ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(0)->get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(0)->get_level(), pl_edge);
      TEST_CHECK_EQUAL(p0.comm_halos.at(0)->get_other(), 1ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(0)->size(), 1ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(0)->get_element(0), 5ul);

      TEST_CHECK_EQUAL(p0.comm_halos.at(1)->get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(1)->get_level(), pl_vertex);
      TEST_CHECK_EQUAL(p0.comm_halos.at(1)->get_other(), 1ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(1)->size(), 1ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(1)->get_element(0), 5ul);

      TEST_CHECK_EQUAL(p0.comm_halos.at(2)->get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(2)->get_level(), pl_edge);
      TEST_CHECK_EQUAL(p0.comm_halos.at(2)->get_other(), 2ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(2)->size(), 1ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(2)->get_element(0), 6ul);

      TEST_CHECK_EQUAL(p0.comm_halos.at(3)->get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(3)->get_level(), pl_vertex);
      TEST_CHECK_EQUAL(p0.comm_halos.at(3)->get_other(), 2ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(3)->size(), 1ul);
      TEST_CHECK_EQUAL(p0.comm_halos.at(3)->get_element(0), 5ul);

      //(3)boundaries for process 0 (should be three edge Halo<0>)
      TEST_CHECK_EQUAL(p0.boundaries.size(), 3ul);

      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_level(), pl_edge);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).size(), 2ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_element(0), 0ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_element(1), 3ul);

      TEST_CHECK_EQUAL(p0.boundaries.at(1).get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(1).get_level(), pl_edge);
      TEST_CHECK_EQUAL(p0.boundaries.at(1).size(), 1ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(1).get_element(0), 1ul);

      TEST_CHECK_EQUAL(p0.boundaries.at(2).get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(2).get_level(), pl_edge);
      TEST_CHECK_EQUAL(p0.boundaries.at(2).size(), 1ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(2).get_element(0), 2ul);

      //(4)attributes
      TEST_CHECK_EQUAL(p0.attrs.at(0)->size(), 6ul);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(0).get())->at(0), 0.);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(1).get())->at(0), 0.);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(0).get())->at(1), 0.5);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(1).get())->at(1), 0.);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(0).get())->at(2), 0.);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(1).get())->at(2), 0.5);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(0).get())->at(3), 1.);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(1).get())->at(3), 0.);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(0).get())->at(4), 1.);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(1).get())->at(4), 0.5);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(0).get())->at(5), 0.5);
      TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(1).get())->at(5), 0.5);
    }
};
PartitioningTest2D<Mem::Main, Index, Algo::Generic, std::vector, std::vector<Index> > ref_test_cpu_v_v_2d("std::vector, std::vector");
PartitioningTest2D<Mem::Main, Index, Algo::Generic, std::vector, std::deque<Index> > ref_test_cpu_v_d_2d("std::vector, std::deque");
PartitioningTest2D<Mem::Main, Index, Algo::Generic, std::deque, std::vector<Index> > ref_test_cpu_d_v_2d("std::deque, std::vector");
PartitioningTest2D<Mem::Main, Index, Algo::Generic, std::deque, std::deque<Index> > ref_test_cpu_d_d_2d("std::deque, std::deque");
