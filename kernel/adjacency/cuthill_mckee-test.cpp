#include <test_system/test_system.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <kernel/adjacency/cuthill_mckee.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Adjacency;

/**
 * \brief Test class for the Cuthill McKee class.
 *
 * \test Tests the Cuthill McKee class.
 *
 * \author Constantin Christof
 */

class CuthillMcKeeTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  CuthillMcKeeTest() :
    TaggedTest<Archs::None, Archs::None>("cuthill_mckee_test")
  {
  }

  bool test_default_root(Graph& g) const
  {
    // Permutation pointer and permutation-array pointer
    Index* perm_array;

    // create Cuthill-McKee object (descending degrees)
    CuthillMcKee cuth1(g, false, CuthillMcKee::root_default,
                      CuthillMcKee::sort_desc);

    // create Cuthill-McKee object (ascending degrees)
    CuthillMcKee cuth2(g, false, CuthillMcKee::root_default,
                      CuthillMcKee::sort_asc);

    // create Cuthill-McKee object (no sorting)
    CuthillMcKee cuth3(g, false, CuthillMcKee::root_default,
                      CuthillMcKee::sort_default);

    // analytic solution (descending degrees)
    Index ref1[12] =
    {
      0, 5, 2, 9, 10, 7, 1, 8, 11, 6, 3, 4
    };

    // analytic solution (ascending degrees)
    Index ref2[12] =
    {
      0, 2, 5, 1, 7, 10, 9, 6, 11, 8, 3, 4
    };

    // analytic solution (no sorting)
    Index ref3[12] =
    {
      0, 2, 5, 1, 7, 9, 10, 6, 8, 11, 3, 4
    };

    // check the permutation arrays
    perm_array = cuth1.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref1[i])
      {
        return false;
      }
    }

    perm_array = cuth2.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref2[i])
      {
        return false;
      }
    }

    perm_array = cuth3.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref3[i])
      {
        return false;
      }
    }
    return true;
  }

  bool test_min_root(Graph& g) const
  {
    // Permutation pointer and permutation-array pointer
    Index* perm_array;

    // create Cuthill-McKee object (descending degrees)
    CuthillMcKee cuth1(g, false, CuthillMcKee::root_minimum_degree,
                      CuthillMcKee::sort_desc);

    // create Cuthill-McKee object (ascending degrees)
    CuthillMcKee cuth2(g, false, CuthillMcKee::root_minimum_degree,
                      CuthillMcKee::sort_asc);

    // create Cuthill-McKee object (no sorting)
    CuthillMcKee cuth3(g, false, CuthillMcKee::root_minimum_degree,
                      CuthillMcKee::sort_default);

    // analytic solution (descending degrees)
    Index ref1[12] =
    {
      4, 3, 6, 1, 2, 0, 5, 9, 10, 7, 8, 11
    };

    // analytic solution (ascending degrees)
    Index ref2[12] =
    {
      4, 3, 6, 1, 2, 0, 5, 7, 10, 9, 11, 8
    };

    // analytic solution (no sorting)
    Index ref3[12] =
    {
      4, 3, 6, 1, 2, 0, 5, 7, 9, 10, 8, 11
    };

    // check the permutation arrays
    perm_array = cuth1.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref1[i])
      {
        return false;
      }
    }

    perm_array = cuth2.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref2[i])
      {
        return false;
      }
    }

    perm_array = cuth3.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref3[i])
      {
        return false;
      }
    }
    return true;
  }

  bool test_max_root(Graph& g) const
  {
    // Permutation pointer and permutation array pointer
    Index* perm_array;

    // create Cuthill-McKee object (descending degrees)
    CuthillMcKee cuth1(g, false, CuthillMcKee::root_maximum_degree,
                      CuthillMcKee::sort_desc);

    // create Cuthill-McKee object (ascending degrees)
    CuthillMcKee cuth2(g, false, CuthillMcKee::root_maximum_degree,
                      CuthillMcKee::sort_asc);

    // create Cuthill-McKee object (no sorting)
    CuthillMcKee cuth3(g, false, CuthillMcKee::root_maximum_degree,
                      CuthillMcKee::sort_default);

    // analytic solution (descending degrees)
    Index ref1[12] =
    {
      5, 9, 10, 0, 7, 8, 11, 2, 1, 6, 3, 4
    };

    // analytic solution (ascending degrees)
    Index ref2[12] =
    {
      5, 0, 7, 10, 9,2, 11, 8, 1, 6, 3, 4
    };

    // analytic solution (no sorting)
    Index ref3[12] =
    {
      5, 0, 7, 9, 10, 2, 8, 11, 1, 6, 3, 4
    };

    // check the permutation arrays
    perm_array = cuth1.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref1[i])
      {
        return false;
      }
    }

    perm_array = cuth2.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref2[i])
      {
        return false;
      }
    }

    perm_array = cuth3.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref3[i])
      {
        return false;
      }
    }
    return true;
  }

  bool test_reverse(Graph& g) const
  {
    // Permutation pointer and permutation array pointer
    Index* perm_array;

    // create Cuthill-McKee object (descending degrees, default root, reverse)
    CuthillMcKee cuth1(g, true, CuthillMcKee::root_default,
                      CuthillMcKee::sort_desc);

    // create Cuthill-McKee object (ascending degrees, min-root, reverse)
    CuthillMcKee cuth2(g, true, CuthillMcKee::root_maximum_degree,
                      CuthillMcKee::sort_asc);

    // create Cuthill-McKee object (no sorting, max-root, reverse)
    CuthillMcKee cuth3(g, true, CuthillMcKee::root_minimum_degree,
                      CuthillMcKee::sort_default);

    // analytic solution (descending degrees)
    Index ref1[12] =
    {
      3, 6, 11, 8, 1, 7, 10, 9, 2, 5, 0, 4
    };

    // analytic solution (ascending degrees)
    Index ref2[12] =
    {
      3, 6, 1, 8, 11, 2, 9, 10, 7, 0, 5, 4
    };

    // analytic solution (no sorting)
    Index ref3[12] =
    {
      4, 11, 8, 10, 9, 7, 5, 0, 2, 1, 6, 3
    };

    // check the permutation arrays
    perm_array = cuth1.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref1[i])
      {
        return false;
      }
    }

    perm_array = cuth2.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref2[i])
      {
        return false;
      }
    }

    perm_array = cuth3.get_permutation().get_perm_pos();

    for(Index i(0); i<12; ++i)
    {
      if(perm_array[i] != ref3[i])
      {
        return false;
      }
    }
    return true;
  }

  virtual void run() const
  {
    // create an adjacency graph
    //       0  1  2  3  4  5  6  7  8  9 10 11   degree
    //    +---------------------------------------------
    //  0 |  0  .  1  .  .  2  .  .  .  .  .  .  |  3
    //  1 |  .  3  4  .  .  .  5  .  .  .  .  .  |  3
    //  2 |  6  7  8  .  .  .  .  .  .  .  .  .  |  3
    //  3 |  .  .  .  9  .  . 10  .  .  .  .  .  |  2
    //  4 |  .  .  .  . 11  .  .  .  .  .  .  .  |  1
    //  5 | 12  .  .  .  . 13  . 14  . 15 16  .  |  5
    //  6 |  . 17  . 18  .  . 19  .  .  .  .  .  |  3
    //  7 |  .  .  .  .  . 20  . 21 22  .  .  .  |  3
    //  8 |  .  .  .  .  .  .  . 23 24 25 26 27  |  5
    //  9 |  .  .  .  .  . 28  .  . 29 30 31 32  |  5
    // 10 |  .  .  .  .  . 33  .  . 34 35 36  .  |  4
    // 11 |  .  .  .  .  .  .  .  . 37 38  . 39  |  3

    Index g_ptr[13] = {0, 3, 6, 9, 11, 12, 17, 20, 23, 28, 33, 37, 40};
    Index g_idx[40] =
    {
      0, 2, 5,
      1, 2, 6,
      0, 1, 2,
      3, 6,
      4,
      0, 5, 7, 9, 10,
      1, 3, 6,
      5, 7, 8,
      7, 8, 9, 10, 11,
      5, 8, 9, 10, 11,
      5, 8, 9, 10,
      8, 9, 11
    };
    Graph g(12, 12, 40, g_ptr, nullptr, g_idx);

    // test Cuthill McKee algorithm
    TEST_CHECK(test_min_root(g));
    TEST_CHECK(test_max_root(g));
    TEST_CHECK(test_default_root(g));
    TEST_CHECK(test_reverse(g));
  }
} cuthill_mckee_test;