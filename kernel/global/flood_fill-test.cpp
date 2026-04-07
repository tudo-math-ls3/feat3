// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/adjacency/graph.hpp>
#include <kernel/base_header.hpp>
#include <kernel/global/flood_fill.hpp>
#include <kernel/util/dist.hpp>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Global;

class FloodFillTest : public UnitTest
{
public:
  FloodFillTest() : UnitTest("FloodFillTest")
  {
  }

  void run() const override
  {

    Dist::Comm comm(Dist::Comm::world());

    if(comm.size() != 2)
    {
      return;
    }

    // 00--01--04--10--06--07
    // |    |          |    |
    // 02--03--05--11--08--09
    if(comm.rank() == 0)
    {
      const Index f_ptr[7] = {0, 2, 5, 7, 10, 12, 14};
      const Index f_idx[14] = {
        1,
        2, // 0
        0,
        3,
        4, // 1
        0,
        3, // 2
        1,
        2,
        5,  // 3
        1,
        10, // 4
        3,
        11  // 5
      };

      Adjacency::Graph graph(6, 12, 14, f_ptr, f_idx);

      std::vector<Index> first_of_rank{0, 6};
      Adjacency::Graph coloring = FloodFill::run(first_of_rank, graph, comm);

      TEST_CHECK_EQUAL(coloring.get_num_nodes_domain(), 1);
    }
    else
    {
      const Index f_ptr[7] = {0, 3, 5, 8, 10, 12, 14};
      const Index f_idx[14] = {
        7,
        8,
        10, // 6
        6,
        9,  // 7
        6,
        9,
        11, // 8
        7,
        8,  // 9
        4,
        6,  // 10
        5,
        8,  // 11
      };

      Adjacency::Graph graph(6, 12, 14, f_ptr, f_idx);
      std::vector<Index> first_of_rank{0, 6};

      Adjacency::Graph coloring = FloodFill::run(first_of_rank, graph, comm);
      TEST_CHECK_EQUAL(coloring.get_num_nodes_domain(), 1);
    }

    // 00--01--04  10--06--07
    // |    |          |    |
    // 02--03--05  11--08--09
    if(comm.rank() == 0)
    {
      const Index f_ptr[7] = {0, 2, 5, 7, 10, 11, 12};
      const Index f_idx[12] = {
        1,
        2, // 0
        0,
        3,
        4, // 1
        0,
        3, // 2
        1,
        2,
        5, // 3
        1, // 4
        3, // 5
      };

      Adjacency::Graph graph(6, 12, 12, f_ptr, f_idx);

      std::vector<Index> first_of_rank{0, 6};
      Adjacency::Graph coloring = FloodFill::run(first_of_rank, graph, comm);

      TEST_CHECK_EQUAL(coloring.get_num_nodes_domain(), 2);
    }
    else
    {
      const Index f_ptr[7] = {0, 3, 5, 8, 10, 11, 12};
      const Index f_idx[12] = {
        7,
        8,
        10, // 6
        6,
        9,  // 7
        6,
        9,
        11, // 8
        7,
        8,  // 9
        6,  // 10
        8,  // 11
      };

      Adjacency::Graph graph(6, 12, 12, f_ptr, f_idx);
      std::vector<Index> first_of_rank{0, 6};

      Adjacency::Graph coloring = FloodFill::run(first_of_rank, graph, comm);
      TEST_CHECK_EQUAL(coloring.get_num_nodes_domain(), 2);
    }
  }
};

static const FloodFillTest flood_fill_test;
