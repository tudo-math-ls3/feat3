// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <deque>
#include <kernel/adjacency/adjactor.hpp>
#include <kernel/adjacency/base.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/global/flood_fill.hpp>
#include <test_system/test_system.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/parti_scotch.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

#ifdef FEAT_HAVE_SCOTCH

class PartiScotchTest :
  public UnitTest
{
public:
  using MeshType = Geometry::ConformalMesh<Shape::Hypercube<2>>;

  PartiScotchTest() :
    UnitTest("parti_scotch-test")
  {
  }

  void run() const override
  {
    constexpr Index num_parts(4);

    Geometry::RefinedUnitCubeFactory<MeshType> factory(2);
    MeshType mesh(factory);

    const Dist::Comm comm = Dist::Comm::world();
    const auto& verts_at_elem = mesh.template get_index_set<MeshType::shape_dim, 0>();

    // Assign weights
    std::vector<Real> weights(std::size_t(mesh.get_num_elements()), Real(1));
    for(std::size_t i(0); i < weights.size(); ++i)
      weights[i] = Real(1 + (i % 3));


    // Run scotch
    Geometry::PartiScotch partitioner(comm);
    TEST_CHECK(partitioner.execute(verts_at_elem, num_parts, weights));

    const Adjacency::Graph elems_at_part = partitioner.build_elems_at_rank();

    // Partition should assign elements to each rank
    TEST_CHECK_EQUAL(elems_at_part.get_num_nodes_domain(), num_parts);
    // Partition should assign rank to each element
    TEST_CHECK_EQUAL(elems_at_part.get_num_nodes_image(), mesh.get_num_elements());

    Index total_assigned(0);
    std::vector<int> elem_counts(std::size_t(mesh.get_num_elements()), 0);

    for(Index part(0); part < num_parts; ++part)
    {
      auto begin = elems_at_part.image_begin(part);
      auto end = elems_at_part.image_end(part);

      for(auto it = begin; it != end; it++)
      {
        const Index element = *it;

        TEST_CHECK(element < mesh.get_num_elements());
        ++elem_counts[element];
        ++total_assigned;
      }
    }

    // All elements should occur in partition
    TEST_CHECK_EQUAL(total_assigned, mesh.get_num_elements());

    // Each element should be assigned to exactly one rank
    for(int count : elem_counts)
    {
      TEST_CHECK_EQUAL(count, 1);
    }

    // Partitions should be connected

    Adjacency::Graph part_at_element(Adjacency::RenderType::transpose, elems_at_part);

    // Build graph G = (V, E), where V are the mesh elements and (a, b) \in E <=> (a neighbor of b and rank(a) = rank(b))

    std::vector<Index> domain_vec;
    std::vector<Index> image_vec;

    const auto& neighbors = mesh.get_neighbors();
    for(Index element(0); element < mesh.get_num_elements(); element++)
    {
      domain_vec.push_back(image_vec.size());

      Index element_part = *part_at_element.image_begin(element);
      for(int i(0); i < 4; i++)
      {
        Index neighbor = neighbors(element, i);
        Index neighbor_part = *part_at_element.image_begin(neighbor);

        if(element_part == neighbor_part)
        {
          image_vec.push_back(neighbor);
        }
      }
    }
    domain_vec.push_back(image_vec.size());

    Adjacency::Graph g(mesh.get_num_elements(), std::move(domain_vec), std::move(image_vec));

    // Run a flood fill on that graph

    std::vector<Index> ff_coloring(mesh.get_num_elements(), ~Index(0));

    Index color(0);
    while(true)
    {
      // Find uncolored starting node
      Index seed(~Index(0));
      for(Index i(0); i < mesh.get_num_elements(); i++)
      {
        if(ff_coloring[i] == ~Index(0))
        {
          seed = i;
          break;
        }
      }

      if(seed == ~Index(0))
      {
        // No seed found, we are done
        break;
      }

      std::deque<Index> stack;
      stack.push_back(seed);
      ff_coloring[seed] = color;

      while(!stack.empty())
      {
        Index element = stack.back();
        stack.pop_back();

        for(Index neighbor: Adjacency::ImageRange(g, element))
        {
          if(ff_coloring[neighbor] == ~Index(0))
          {
            ff_coloring[neighbor] = color;
            stack.push_back(neighbor);
          }
        }
      }

      color++;
    }

    std::array<std::set<Index>, 4> colors_per_rank;

    for(Index part(0); part < num_parts; part++)
    {
      for(Index element : Adjacency::ImageRange(elems_at_part, part))
      {
        colors_per_rank.at(part).insert(ff_coloring[element]);
      }
    }

    // Check that each partition only contains a single color
    for(Index i(0); i < num_parts; i++)
    {
      TEST_CHECK_EQUAL(colors_per_rank.at(i).size(), 1);
    }
  }
} parti_scotch_test;

#endif
