#include <test_system/test_system.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/colouring.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Adjacency;

/**
 * \brief Test class for the Colouring class.
 *
 * \test Tests the Colouring class.
 *
 * \author Constantin Christof
 */

  class ColouringTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  ColouringTest() :
    TaggedTest<Archs::None, Archs::None>("colouring_test")
  {
  }

  bool test_c(Colouring& c) const
  {
    // check against analytic solution
    Index col_ref[7] = {0, 0, 1, 1, 1, 2, 2};

    // fetch the colouring array
    Index* _colouring = c.get_colouring();

    // compare index arrays
    for(Index j(0); j < 7; ++j)
    {
      if(col_ref[j] != _colouring[j])
      {
        return false;
      }
    }
    return true;
  }

  bool test_c_ordered(Colouring& c) const
  {
    // check against analytic solution
    Index col_ref[7] = {2, 2, 1, 1, 0, 0, 0};

    // fetch the colouring array
    Index* _colouring = c.get_colouring();

    // compare index arrays
    for(Index j(0); j < 7; ++j)
    {
      if(col_ref[j] != _colouring[j])
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
    Graph g(7, 7, 20, g_ptr, nullptr, g_idx);

    // create a colouring object for this graph
    Colouring c(g);

    // validate
    TEST_CHECK(test_c(c));

    // ordered colouring

    // define permutation array
    Index _order[7] =
    {
      5, 2, 4, 1, 3, 0, 6
    };

    // create a colouring object corresponding to _order
    Colouring co(g, _order);

    // validate
    TEST_CHECK(test_c_ordered(co));

  }

} colouring_test;
