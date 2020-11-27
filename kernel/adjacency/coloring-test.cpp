// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/coloring.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Adjacency;

/**
 * \brief Test class for the Coloring class.
 *
 * \test Tests the Coloring class.
 *
 * \author Constantin Christof
 */

class ColoringTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  ColoringTest() :
    TaggedTest<Archs::None, Archs::None>("coloring_test")
  {
  }

  virtual ~ColoringTest()
  {
  }

  bool test_c(Coloring& c) const
  {
    // check against analytic solution
    Index col_ref[7] = {0, 0, 1, 1, 1, 2, 2};

    // fetch the coloring array
    Index* _coloring = c.get_coloring();

    // compare index arrays
    for(Index j(0); j < 7; ++j)
    {
      if(col_ref[j] != _coloring[j])
      {
        return false;
      }
    }
    return true;
  }

  bool test_c_ordered(Coloring& c) const
  {
    // check against analytic solution
    Index col_ref[7] = {2, 2, 1, 1, 0, 0, 0};

    // fetch the coloring array
    Index* _coloring = c.get_coloring();

    // compare index arrays
    for(Index j(0); j < 7; ++j)
    {
      if(col_ref[j] != _coloring[j])
      {
        return false;
      }
    }
    return true;
  }

  virtual void run() const override
  {
    // create an adjacency graph
    //      0  1  2  3  4  5  6
    //   +---------------------
    // 0 |  0  .  1  2  .  3  .
    // 1 |  .  .  4  .  5  .  6
    // 2 |  7  8  9  .  . 10  .
    // 3 | 11  .  . 12  .  . 13
    // 4 |  . 14  .  .  .  .  .
    // 5 | 15  . 16  .  .  .  .
    // 6 |  . 17  . 18  .  . 19

    Index g_ptr[8] = {0, 4, 7, 11, 14, 15, 17, 20};
    Index g_idx[20] =
    {
      0, 2, 3, 5,
      2, 4, 6,
      0, 1, 2, 5,
      0, 3, 6,
      1,
      0, 2,
      1, 3, 6
    };
    Graph g(7, 7, 20, g_ptr, g_idx);

    // create a coloring object for this graph
    Coloring c(g);

    // validate
    TEST_CHECK(test_c(c));

    // ordered coloring

    // define permutation array
    Index _order[7] =
    {
      5, 2, 4, 1, 3, 0, 6
    };

    // create a coloring object corresponding to _order
    Coloring co(g, _order);

    // validate
    TEST_CHECK(test_c_ordered(co));

  }
} coloring_test;
