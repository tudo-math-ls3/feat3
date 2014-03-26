#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/partitioning.hpp>
#include<kernel/foundation/aura.hpp>
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

      OT_<std::shared_ptr<HaloBase<Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> > > > halos;
      OT_<Halo<0, PLVertex, Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>, std::allocator<Halo<0, PLVertex, Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> > > boundaries;

      boundaries.push_back(Halo<0, PLVertex, Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLVertex, Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
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

      OT_<Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>, std::allocator<Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> > > boundaries;
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
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
      TEST_CHECK_EQUAL(sm_edge_6_adj_vert.at(0), 5ul);
      TEST_CHECK_EQUAL(sm_edge_6_adj_vert.at(1), 4ul);

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

template<typename Tag_, typename IndexType_, typename Algo_, template<typename, typename> class OT_, typename IT_>
class PartitioningTest3D:
  public TaggedTest<Tag_, IndexType_, Algo_>
{
  public:
    PartitioningTest3D(const std::string & tag) :
      TaggedTest<Tag_, Index, Algo_>("PartitioningTest3D<" + tag + ">")
    {
    }

    virtual void run() const
    {
      /*
       * (0,1,1)  (1,1,1)
       *      *----*
       *     /    /|
       *(0,1,0)(1,1,0)
       *   *----*  *(1,0,1)
       *   | /  | /
       *   |/   |/
       *   *----*
       *(0,0,0) (1,0,0)
       */

      //create attributes for vertex coords
      OT_<Attribute<double, OT_>, std::allocator<Attribute<double, OT_> > > attrs;
      attrs.push_back(Attribute<double, OT_>()); //vertex x-coords
      attrs.push_back(Attribute<double, OT_>()); //vertex y-coords
      attrs.push_back(Attribute<double, OT_>()); //vertex z-coords

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(0));
      attrs.at(2).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(0));
      attrs.at(2).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(1));
      attrs.at(2).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(1));
      attrs.at(2).get_data().push_back(double(0));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(0));
      attrs.at(2).get_data().push_back(double(1));

      attrs.at(0).get_data().push_back(double(1));
      attrs.at(1).get_data().push_back(double(1));
      attrs.at(2).get_data().push_back(double(1));

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(0));
      attrs.at(2).get_data().push_back(double(1));

      attrs.at(0).get_data().push_back(double(0));
      attrs.at(1).get_data().push_back(double(1));
      attrs.at(2).get_data().push_back(double(1));

      /*  2    3
       *  *-1--*
       *  2    |
       *  |    3
       *  *--0-*
       *  0    1
       */

      //creating foundation mesh
      Foundation::Mesh<Dim3D, Foundation::Topology<IndexType_, OT_, IT_>, OT_> m(0);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);
      m.add_polytope(pl_vertex);

      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);
      m.add_polytope(pl_edge);

      m.add_polytope(pl_face);
      m.add_polytope(pl_face);
      m.add_polytope(pl_face);
      m.add_polytope(pl_face);
      m.add_polytope(pl_face);
      m.add_polytope(pl_face);

      m.add_polytope(pl_polyhedron);

      //vertex 0
      m.add_adjacency(pl_vertex, pl_edge, 0, 0);
      m.add_adjacency(pl_vertex, pl_edge, 0, 2);
      m.add_adjacency(pl_vertex, pl_edge, 0, 10);

      m.add_adjacency(pl_vertex, pl_face, 0, 0);
      m.add_adjacency(pl_vertex, pl_face, 0, 3);
      m.add_adjacency(pl_vertex, pl_face, 0, 4);

      m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

      //vertex 1
      m.add_adjacency(pl_vertex, pl_edge, 1, 0);
      m.add_adjacency(pl_vertex, pl_edge, 1, 3);
      m.add_adjacency(pl_vertex, pl_edge, 1, 4);

      m.add_adjacency(pl_vertex, pl_face, 1, 0);
      m.add_adjacency(pl_vertex, pl_face, 1, 2);
      m.add_adjacency(pl_vertex, pl_face, 1, 4);

      m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

      //vertex 2
      m.add_adjacency(pl_vertex, pl_edge, 2, 1);
      m.add_adjacency(pl_vertex, pl_edge, 2, 2);
      m.add_adjacency(pl_vertex, pl_edge, 2, 11);

      m.add_adjacency(pl_vertex, pl_face, 2, 0);
      m.add_adjacency(pl_vertex, pl_face, 2, 3);
      m.add_adjacency(pl_vertex, pl_face, 2, 5);

      m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

      //vertex 3
      m.add_adjacency(pl_vertex, pl_edge, 3, 1);
      m.add_adjacency(pl_vertex, pl_edge, 3, 3);
      m.add_adjacency(pl_vertex, pl_edge, 3, 5);

      m.add_adjacency(pl_vertex, pl_face, 3, 2);
      m.add_adjacency(pl_vertex, pl_face, 3, 0);
      m.add_adjacency(pl_vertex, pl_face, 3, 5);

      m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

      //vertex 4
      m.add_adjacency(pl_vertex, pl_edge, 4, 4);
      m.add_adjacency(pl_vertex, pl_edge, 4, 6);
      m.add_adjacency(pl_vertex, pl_edge, 4, 7);

      m.add_adjacency(pl_vertex, pl_face, 4, 1);
      m.add_adjacency(pl_vertex, pl_face, 4, 2);
      m.add_adjacency(pl_vertex, pl_face, 4, 4);

      m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

      //vertex 5
      m.add_adjacency(pl_vertex, pl_edge, 5, 5);
      m.add_adjacency(pl_vertex, pl_edge, 5, 6);
      m.add_adjacency(pl_vertex, pl_edge, 5, 8);

      m.add_adjacency(pl_vertex, pl_face, 5, 1);
      m.add_adjacency(pl_vertex, pl_face, 5, 2);
      m.add_adjacency(pl_vertex, pl_face, 5, 5);

      m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

      //vertex 6
      m.add_adjacency(pl_vertex, pl_edge, 6, 7);
      m.add_adjacency(pl_vertex, pl_edge, 6, 9);
      m.add_adjacency(pl_vertex, pl_edge, 6, 10);

      m.add_adjacency(pl_vertex, pl_face, 6, 1);
      m.add_adjacency(pl_vertex, pl_face, 6, 3);
      m.add_adjacency(pl_vertex, pl_face, 6, 4);

      m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

      //vertex 7
      m.add_adjacency(pl_vertex, pl_edge, 7, 8);
      m.add_adjacency(pl_vertex, pl_edge, 7, 9);
      m.add_adjacency(pl_vertex, pl_edge, 7, 11);

      m.add_adjacency(pl_vertex, pl_face, 7, 1);
      m.add_adjacency(pl_vertex, pl_face, 7, 3);
      m.add_adjacency(pl_vertex, pl_face, 7, 5);

      m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

      Foundation::Mesh<Dim3D, Foundation::Topology<IndexType_, OT_, IT_>, OT_> m_fine(m);

      OT_<Halo<0, PLFace, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>, std::allocator<Halo<0, PLFace, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> > > boundaries;
      boundaries.push_back(Halo<0, PLFace, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLFace, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLFace, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLFace, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLFace, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLFace, Mesh<Dim3D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.at(0).push_back(0);
      boundaries.at(1).push_back(1);
      boundaries.at(2).push_back(2);
      boundaries.at(3).push_back(3);
      boundaries.at(3).push_back(4);
      boundaries.at(3).push_back(5);

      Index num_procs(3);
      Index rank(0);
      //Index level(4);

      PData<Dim3D, Topology<IndexType_, OT_, IT_>, OT_, Mesh, double> p0(Partitioning<Tag_,
                                                                                      Algo_,
                                                                                      Dim3D,
                                                                                      0,
                                                                                      pl_vertex>::execute(m,
                                                                                                          boundaries,
                                                                                                          num_procs, rank,
                                                                                                          attrs
                                                                                                          ));

      //inspect p0 for (1)submesh data, (2)comm_halos, (3)boundaries, (4)attributes
      //(1)submesh data on process 0
      TEST_CHECK_EQUAL(p0.submesh->get_topologies().at(ipi_polyhedron_vertex).size(), 3ul);
      TEST_CHECK_EQUAL(p0.submesh->get_topologies().at(ipi_face_vertex).size(), 16ul);
      TEST_CHECK_EQUAL(p0.submesh->get_topologies().at(ipi_edge_vertex).size(), 28ul);
      TEST_CHECK_EQUAL(p0.submesh->get_topologies().at(ipi_vertex_edge).size(), 16ul);

      TEST_CHECK_EQUAL(p0.submesh->get_map().size(), 3ul);
      TEST_CHECK_EQUAL(p0.submesh->get_map().at(0), 0ul);
      TEST_CHECK_EQUAL(p0.submesh->get_map().at(1), 1ul);

      auto sm_polyhedron_0_adj_vert(p0.submesh->get_adjacent_polytopes(pl_polyhedron, pl_vertex, 0));
      TEST_CHECK_EQUAL(sm_polyhedron_0_adj_vert.size(), 8ul);
      TEST_CHECK_EQUAL(sm_polyhedron_0_adj_vert.at(0), 0ul);
      TEST_CHECK_EQUAL(sm_polyhedron_0_adj_vert.at(1), 1ul);
      TEST_CHECK_EQUAL(sm_polyhedron_0_adj_vert.at(2), 2ul);
      TEST_CHECK_EQUAL(sm_polyhedron_0_adj_vert.at(3), 8ul);
      TEST_CHECK_EQUAL(sm_polyhedron_0_adj_vert.at(4), 10ul);
      TEST_CHECK_EQUAL(sm_polyhedron_0_adj_vert.at(5), 12ul);
      TEST_CHECK_EQUAL(sm_polyhedron_0_adj_vert.at(6), 13ul);
      TEST_CHECK_EQUAL(sm_polyhedron_0_adj_vert.at(7), 15ul);

      auto sm_polyhedron_1_adj_vert(p0.submesh->get_adjacent_polytopes(pl_polyhedron, pl_vertex, 1));
      TEST_CHECK_EQUAL(sm_polyhedron_1_adj_vert.size(), 8ul);
      TEST_CHECK_EQUAL(sm_polyhedron_1_adj_vert.at(0), 3ul);
      TEST_CHECK_EQUAL(sm_polyhedron_1_adj_vert.at(1), 1ul);
      TEST_CHECK_EQUAL(sm_polyhedron_1_adj_vert.at(2), 4ul);
      TEST_CHECK_EQUAL(sm_polyhedron_1_adj_vert.at(3), 5ul);
      TEST_CHECK_EQUAL(sm_polyhedron_1_adj_vert.at(4), 10ul);
      TEST_CHECK_EQUAL(sm_polyhedron_1_adj_vert.at(5), 11ul);
      TEST_CHECK_EQUAL(sm_polyhedron_1_adj_vert.at(6), 13ul);
      TEST_CHECK_EQUAL(sm_polyhedron_1_adj_vert.at(7), 15ul);

      auto sm_polyhedron_2_adj_vert(p0.submesh->get_adjacent_polytopes(pl_polyhedron, pl_vertex, 2));
      TEST_CHECK_EQUAL(sm_polyhedron_2_adj_vert.size(), 8ul);
      TEST_CHECK_EQUAL(sm_polyhedron_2_adj_vert.at(0), 6ul);
      TEST_CHECK_EQUAL(sm_polyhedron_2_adj_vert.at(1), 9ul);
      TEST_CHECK_EQUAL(sm_polyhedron_2_adj_vert.at(2), 4ul);
      TEST_CHECK_EQUAL(sm_polyhedron_2_adj_vert.at(3), 7ul);
      TEST_CHECK_EQUAL(sm_polyhedron_2_adj_vert.at(4), 10ul);
      TEST_CHECK_EQUAL(sm_polyhedron_2_adj_vert.at(5), 11ul);
      TEST_CHECK_EQUAL(sm_polyhedron_2_adj_vert.at(6), 14ul);
      TEST_CHECK_EQUAL(sm_polyhedron_2_adj_vert.at(7), 15ul);

      //intersection in {1,10,13,15} implies comm_intersection on vertex level must be {1,10,13,15}
      auto vertex_comm_int_for_polyhedra_0_1(p0.submesh->get_comm_intersection(pl_polyhedron, pl_vertex, 0, 1));
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_0_1.size(), 4ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_0_1.at(0), 1ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_0_1.at(1), 10ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_0_1.at(2), 13ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_0_1.at(3), 15ul);
      //intersection in {10,15} implies comm_intersection on vertex level must be {10,15}
      auto vertex_comm_int_for_polyhedra_0_2(p0.submesh->get_comm_intersection(pl_polyhedron, pl_vertex, 0, 2));
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_0_2.size(), 2ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_0_2.at(0), 10ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_0_2.at(1), 15ul);
      //intersection in {4,10,11,15} implies comm_intersection on vertex level must be {4,10,11,15}
      auto vertex_comm_int_for_polyhedra_1_2(p0.submesh->get_comm_intersection(pl_polyhedron, pl_vertex, 1, 2));
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_1_2.size(), 4ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_1_2.at(0), 4ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_1_2.at(1), 10ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_1_2.at(2), 11ul);
      TEST_CHECK_EQUAL(vertex_comm_int_for_polyhedra_1_2.at(3), 15ul);

      const Index ref_face_adj_vert[16][4] =
      {
        { 0, 1, 2, 10 },
        { 3, 4, 5, 11},
        { 0, 2, 8, 12 },
        { 0, 1, 8, 13 },
        { 1, 3, 10, 4 },
        { 10, 4, 9, 6 },
        { 4, 6, 11, 7 },
        { 1, 3, 13, 5},
        { 9, 6, 14, 7 },
        { 15, 10, 11, 4 },
        { 15, 10, 12, 2 },
        { 15, 10, 13, 1 },
        { 15, 10, 14, 9 },
        { 15, 11, 13, 5 },
        { 15, 11, 14, 7 },
        { 15, 12, 13, 8 }
      };
      for( Index j(0); j < 16; j++)
      {
        auto sm_face_j_adj_vert(p0.submesh->get_adjacent_polytopes(pl_face, pl_vertex, j));
        TEST_CHECK_EQUAL(sm_face_j_adj_vert.size(), 4ul);
        TEST_CHECK_EQUAL(sm_face_j_adj_vert.at(0), ref_face_adj_vert[j][0]);
        TEST_CHECK_EQUAL(sm_face_j_adj_vert.at(1), ref_face_adj_vert[j][1]);
        TEST_CHECK_EQUAL(sm_face_j_adj_vert.at(2), ref_face_adj_vert[j][2]);
        TEST_CHECK_EQUAL(sm_face_j_adj_vert.at(3), ref_face_adj_vert[j][3]);
      }

// => comm_intersection on face level
      auto face_comm_int_for_polyhedra_0_1(p0.submesh->get_comm_intersection(pl_polyhedron, pl_face, 0, 1));
      TEST_CHECK_EQUAL(face_comm_int_for_polyhedra_0_1.size(), 1ul);
      TEST_CHECK_EQUAL(face_comm_int_for_polyhedra_0_1.at(0), 11ul);
      auto face_comm_int_for_polyhedra_0_2(p0.submesh->get_comm_intersection(pl_polyhedron, pl_face, 0, 2));
      TEST_CHECK_EQUAL(face_comm_int_for_polyhedra_0_2.size(), 0ul);
      auto face_comm_int_for_polyhedra_1_2(p0.submesh->get_comm_intersection(pl_polyhedron, pl_face, 1, 2));
      TEST_CHECK_EQUAL(face_comm_int_for_polyhedra_1_2.size(), 1ul);
      TEST_CHECK_EQUAL(face_comm_int_for_polyhedra_1_2.at(0), 9ul);


      const Index ref_edge_adj_vert[28][2] =
      {
        { 0, 1 },
        { 0, 2 },
        { 3, 4 },  { 3, 5 },
        { 6, 7 },
        { 0, 8 },
        { 1, 3 },
        { 9, 6 },
        { 4, 6 },
        { 1, 10 },
        { 9, 10 },
        { 2, 10 },
        { 4, 10 },
        { 4, 11 },
        { 5, 11 },
        { 7, 11 },
        { 2, 12 },
        { 8, 12 },
        { 1, 13},
        { 5, 13 },
        { 8, 13 },
        { 9, 14 },
        { 7, 14 },
        { 10, 15 },
        { 11, 15 },
        { 12, 15 },
        { 13, 15 },
        { 14, 15 }
      };
      for (Index j(0); j < p0.submesh->get_topologies().at(ipi_edge_vertex).size(); j++)
      {
        auto sm_edge_j_adj_vert(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, j));
      TEST_CHECK_EQUAL(sm_edge_j_adj_vert.size(), 2ul);
      TEST_CHECK_EQUAL(sm_edge_j_adj_vert.at(0), ref_edge_adj_vert[j][0]);
      TEST_CHECK_EQUAL(sm_edge_j_adj_vert.at(1), ref_edge_adj_vert[j][1]);
      }


      //=> comm_intersection on edge level
      auto edge_comm_int_for_polyhedra_0_1(p0.submesh->get_comm_intersection(pl_polyhedron, pl_edge, 0, 1));
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_0_1.size(), 4ul);
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_0_1.at(0), 9ul);
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_0_1.at(1), 18ul);
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_0_1.at(2), 23ul);
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_0_1.at(3), 26ul);
      auto edge_comm_int_for_polyhedra_0_2(p0.submesh->get_comm_intersection(pl_polyhedron, pl_edge, 0, 2));
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_0_2.size(), 1ul);
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_0_2.at(0), 23ul);
      auto edge_comm_int_for_polyhedra_1_2(p0.submesh->get_comm_intersection(pl_polyhedron, pl_edge, 1, 2));
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_1_2.size(), 4ul);
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_1_2.at(0), 12ul);
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_1_2.at(1), 13ul);
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_1_2.at(2), 23ul);
      TEST_CHECK_EQUAL(edge_comm_int_for_polyhedra_1_2.at(3), 24ul);

      //(2)comm halos for process 0
      const Index ref_level[15] = {2, 2, 2, 1, 1, 1, 1, 0, 0, 2, 2, 1, 1, 1, 0};
      const Index ref_other[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2};
      const Index ref_element[15] = {15, 13, 14, 26, 24, 26, 24, 15, 15, 10, 12, 25, 23, 27, 15};
      for( Index j(0) ; j < 15 ; j++)
      {
        TEST_CHECK_EQUAL(p0.comm_halos.at(j)->get_overlap(), 0ul);
        TEST_CHECK_EQUAL(Index(p0.comm_halos.at(j)->get_level()), ref_level[j]);
        TEST_CHECK_EQUAL(p0.comm_halos.at(j)->get_other(), ref_other[j]);
        TEST_CHECK_EQUAL(p0.comm_halos.at(j)->size(), 1ul);
        TEST_CHECK_EQUAL(p0.comm_halos.at(j)->get_element(0), ref_element[j]);
      }

      //(3)boundaries for process 0 (should be three edge Halo<0>)
      TEST_CHECK_EQUAL(p0.boundaries.size(), 3ul);

      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_level(), pl_face);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).size(), 3ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_element(0), 0ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_element(1), 5ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(0).get_element(2), 4ul);

      TEST_CHECK_EQUAL(p0.boundaries.at(1).get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(1).get_level(), pl_face);
      TEST_CHECK_EQUAL(p0.boundaries.at(1).size(), 2ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(1).get_element(0), 1ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(1).get_element(1), 6ul);

      TEST_CHECK_EQUAL(p0.boundaries.at(2).get_overlap(), 0ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(2).get_level(), pl_face);
      TEST_CHECK_EQUAL(p0.boundaries.at(2).size(), 4ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(2).get_element(0), 2ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(2).get_element(1), 3ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(2).get_element(2), 7ul);
      TEST_CHECK_EQUAL(p0.boundaries.at(2).get_element(3), 8ul);

      //(4)attributes
      const double ref_attrs[16][3] =
      {
        {0., 0., 0.},
        {0.5, 0., 0.},
        {0., 0.5, 0.},
        {1., 0., 0.},
        {1., 0.5, 0.},
        {1., 0., 0.5},
        {1., 1., 0.},
        {1., 1., 0.5},
        {0., 0., 0.5},
        {0.5, 1., 0.},
        {0.5, 0.5, 0.},
        {1., 0.5, 0.5},
        {0., 0.5, 0.5},
        {0.5, 0., 0.5},
        {0.5, 1., 0.5},
        {0.5, 0.5, 0.5}
      };
      TEST_CHECK_EQUAL(p0.attrs.at(0)->size(), 16ul);
      for( Index j(0); j < 16; j++ )
      {
        TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(0).get())->at(j), ref_attrs[j][0]);
        TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(1).get())->at(j), ref_attrs[j][1]);
        TEST_CHECK_EQUAL(((Foundation::Attribute<double, OT_>*)p0.attrs.at(2).get())->at(j), ref_attrs[j][2]);
      }

    }
};
PartitioningTest3D<Mem::Main, Index, Algo::Generic, std::vector, std::vector<Index> > ref_test_cpu_v_v_3d("std::vector, std::vector");
PartitioningTest3D<Mem::Main, Index, Algo::Generic, std::vector, std::deque<Index> > ref_test_cpu_v_d_3d("std::vector, std::deque");
PartitioningTest3D<Mem::Main, Index, Algo::Generic, std::deque, std::vector<Index> > ref_test_cpu_d_v_3d("std::deque, std::vector");
PartitioningTest3D<Mem::Main, Index, Algo::Generic, std::deque, std::deque<Index> > ref_test_cpu_d_d_3d("std::deque, std::deque");

template<typename Tag_, typename IndexType_, typename Algo_, template<typename, typename> class OT_, typename IT_>
class PartitioningTest2DHaloOrdering:
  public TaggedTest<Tag_, IndexType_, Algo_>
{
  public:
    PartitioningTest2DHaloOrdering(const std::string & tag) :
      TaggedTest<Tag_, Index, Algo_>("PartitioningTest2DHaloOrdering<" + tag + ">")
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

      OT_<Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>, std::allocator<Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> > > boundaries;
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m));
      boundaries.at(0).push_back(0);
      boundaries.at(1).push_back(1);
      boundaries.at(2).push_back(2);
      boundaries.at(3).push_back(3);

      PData<Dim2D, Topology<IndexType_, OT_, IT_>, OT_, Mesh, double> p0(Partitioning<Tag_,
                                                                                      Algo_,
                                                                                      Dim2D,
                                                                                      0,
                                                                                      pl_vertex>::execute(m,
                                                                                                          boundaries,
                                                                                                          2, 0,
                                                                                                          attrs
                                                                                                          ));
      PData<Dim2D, Topology<IndexType_, OT_, IT_>, OT_, Mesh, double> p1(Partitioning<Tag_,
                                                                                      Algo_,
                                                                                      Dim2D,
                                                                                      0,
                                                                                      pl_vertex>::execute(m,
                                                                                                          boundaries,
                                                                                                          2, 1,
                                                                                                          attrs
                                                                                                          ));

      for(auto& h0_i : p0.comm_halos)
        for(auto& h1_i : p1.comm_halos)
          if(h0_i->get_other() == h1_i->get_other())
          {
            //every halo element must appear with same index in both patches' submeshes' maps
            for(Index i(0) ; i < h0_i->size() ; ++i)
            {
              Index elem_nr_loc_0(h0_i->get_element(i));
              Index elem_nr_loc_1(h1_i->get_element(i));

              auto adj_v_0(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, elem_nr_loc_0));
              auto adj_v_1(p1.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, elem_nr_loc_1));
              TEST_CHECK_EQUAL(
                  ((Attribute<double, OT_>*)(p0.attrs.at(0).get()))->get_data().at(adj_v_0.at(0)),
                  ((Attribute<double, OT_>*)(p1.attrs.at(0).get()))->get_data().at(adj_v_1.at(0))
                  ); //x coords of first vertex
              TEST_CHECK_EQUAL(
                  ((Attribute<double, OT_>*)(p0.attrs.at(0).get()))->get_data().at(adj_v_0.at(1)),
                  ((Attribute<double, OT_>*)(p1.attrs.at(0).get()))->get_data().at(adj_v_1.at(1))
                  ); //x coords of second vertex
              TEST_CHECK_EQUAL(
                  ((Attribute<double, OT_>*)(p0.attrs.at(1).get()))->get_data().at(adj_v_0.at(0)),
                  ((Attribute<double, OT_>*)(p1.attrs.at(1).get()))->get_data().at(adj_v_1.at(0))
                  ); //y coords of first vertex
              TEST_CHECK_EQUAL(
                  ((Attribute<double, OT_>*)(p0.attrs.at(1).get()))->get_data().at(adj_v_0.at(1)),
                  ((Attribute<double, OT_>*)(p1.attrs.at(1).get()))->get_data().at(adj_v_1.at(1))
                  ); //y coords of second vertex
            }
          }

      auto h_new0(Aura<Mem::Main, Algo::Generic, Halo<0, PLEdge, Mesh<Dim2D> > >::value(p0.comm_halos));
      auto h_new1(Aura<Mem::Main, Algo::Generic, Halo<0, PLEdge, Mesh<Dim2D> > >::value(p1.comm_halos));
      for(Index i(0) ; i < h_new0.size() ; ++i)
      {
        Index elem_nr_loc_0(h_new0.get_element(i));
        Index elem_nr_loc_1(h_new1.get_element(i));

        auto adj_v_0(p0.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, elem_nr_loc_0));
        auto adj_v_1(p1.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, elem_nr_loc_1));
        TEST_CHECK_EQUAL(
            ((Attribute<double, OT_>*)(p0.attrs.at(0).get()))->get_data().at(adj_v_0.at(0)),
            ((Attribute<double, OT_>*)(p1.attrs.at(0).get()))->get_data().at(adj_v_1.at(0))
            ); //x coords of first vertex
        TEST_CHECK_EQUAL(
            ((Attribute<double, OT_>*)(p0.attrs.at(0).get()))->get_data().at(adj_v_0.at(1)),
            ((Attribute<double, OT_>*)(p1.attrs.at(0).get()))->get_data().at(adj_v_1.at(1))
            ); //x coords of second vertex
        TEST_CHECK_EQUAL(
            ((Attribute<double, OT_>*)(p0.attrs.at(1).get()))->get_data().at(adj_v_0.at(0)),
            ((Attribute<double, OT_>*)(p1.attrs.at(1).get()))->get_data().at(adj_v_1.at(0))
            ); //y coords of first vertex
        TEST_CHECK_EQUAL(
            ((Attribute<double, OT_>*)(p0.attrs.at(1).get()))->get_data().at(adj_v_0.at(1)),
            ((Attribute<double, OT_>*)(p1.attrs.at(1).get()))->get_data().at(adj_v_1.at(1))
            ); //y coords of second vertex
      }
    }
};
PartitioningTest2DHaloOrdering<Mem::Main, Index, Algo::Generic, std::vector, std::vector<Index> > parthaloorder_test_cpu_v_v_2d("std::vector, std::vector");
