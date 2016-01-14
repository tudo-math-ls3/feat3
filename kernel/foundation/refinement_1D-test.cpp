#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/refinement.hpp>
#include<kernel/foundation/export.hpp>
#include<kernel/archs.hpp>

#include<kernel/foundation/mesh_control.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/export_vtk.hpp>

#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class RefinementTest1D:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    RefinementTest1D(const std::string & tag) :
      TaggedTest<Tag_, Index>("RefinementTest1D<" + tag + ">")
    {
    }

    virtual void run() const override
    {
      /* (0)  (1)
       *  *----*
       */

      //create attributes for vertex coords
      OT_< Attribute<double, OT_>, std::allocator<Attribute<double, OT_> > > attrs;
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


      Foundation::Mesh<Dim1D, Foundation::Topology<IndexType_, OT_, IT_>, OT_> m_fine(m);

      //set up halos
      OT_<std::shared_ptr<HaloBase<Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> > > > halos;
      halos.push_back(std::shared_ptr<HaloBase<Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_> >(new Halo<1, PLEdge, Mesh<Dim1D, Topology<IndexType_, OT_, IT_>, OT_>, double, OT_>(m)));
      halos.at(0)->push_back(0);


      Refinement<Mem::Main,
                 mrt_standard,
                 hrt_refine>::execute(m_fine, &halos, attrs);


      //testing adjacencies 'edges at vertices'
      TEST_CHECK_EQUAL(m_fine.get_topologies().at(ipi_vertex_edge).size(), 3ul);
      TEST_CHECK_EQUAL(m_fine.get_topologies().at(ipi_edge_vertex).size(), 2ul);

      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_vertex, pl_edge, 0).size(), 1ul);
      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_vertex, pl_edge, 0).at(0), 0ul);

      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_vertex, pl_edge, 1).size(), 1ul);
      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_vertex, pl_edge, 1).at(0), 1ul);

      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_vertex, pl_edge, 2).size(), 2ul);
      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_vertex, pl_edge, 2).at(0), 0ul);
      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_vertex, pl_edge, 2).at(1), 1ul);

      //testing adjacencies 'vertices at edges' -> ORIENTATION PROOFNESS
      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 0).at(0), 0ul);
      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 0).at(1), 2ul);

      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 1).at(0), 2ul);
      TEST_CHECK_EQUAL(m_fine.get_adjacent_polytopes(pl_edge, pl_vertex, 1).at(1), 1ul);

      //testing vertex coords
      TEST_CHECK_EQUAL(attrs.at(0).at(0), double(0));
      TEST_CHECK_EQUAL(attrs.at(0).at(1), double(1));
      TEST_CHECK_EQUAL(attrs.at(0).at(2), double(0.5));

      TEST_CHECK_EQUAL(halos.at(0)->get_elements().size(), IndexType_(2));
      TEST_CHECK_EQUAL(halos.at(0)->get_elements().at(0), IndexType_(0));
      TEST_CHECK_EQUAL(halos.at(0)->get_elements().at(1), IndexType_(1));
    }
};
RefinementTest1D<Mem::Main, Index, std::vector, std::vector<Index> > ref_test_cpu_v_v("std::vector, std::vector");
RefinementTest1D<Mem::Main, Index, std::vector, std::deque<Index> > ref_test_cpu_v_d("std::vector, std::deque");
RefinementTest1D<Mem::Main, Index, std::deque, std::vector<Index> > ref_test_cpu_d_v("std::deque, std::vector");
RefinementTest1D<Mem::Main, Index, std::deque, std::deque<Index> > ref_test_cpu_d_d("std::deque, std::deque");
