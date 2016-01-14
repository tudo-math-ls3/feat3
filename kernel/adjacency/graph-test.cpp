#include <test_system/test_system.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/colouring.hpp>
#include <kernel/adjacency/permutation.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Adjacency;

typedef CompositeAdjactor<Graph,Graph> CAGG;

/**
 * \brief Test class for the Graph class.
 *
 * \test Tests the Graph class.
 *
 * \author Peter Zajac
 */
class GraphTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  GraphTest() :
    TaggedTest<Archs::None, Archs::None>("graph_test")
  {
  }

  bool test_f(Graph& f) const
  {
    CONTEXT("GraphTest::test_f()");

    // fetch the graph's arrays
    Index* f_pp = f.get_domain_ptr();
    Index* f_di = f.get_image_idx();

    // check against analytic solution
    //      0  1  2  3  4
    //   +---------------
    // 0 |  0  .  1  .  2
    // 1 |  .  .  3  4  .
    // 2 |  5  6  7  .  .
    // 3 |  8  .  .  9 10
    // 4 |  . 11  .  . 12
    // 5 | 13  . 14  .  .
    // 6 |  . 15  . 16  .

    Index f_ptr[8] = {0, 3, 5, 8, 11, 13, 15, 17};
    Index f_idx[17] =
    {
      0, 2, 4,
      2, 3,
      0, 1, 2,
      0, 3, 4,
      1, 4,
      0, 2,
      1, 3
    };

    // check dimensions
    if(f.get_num_nodes_domain() != 7)
      return false;
    if(f.get_num_nodes_image() != 5)
      return false;

    // check degree
    if(f.degree() != 3)
      return false;

    // compare pointer arrays
    for(int i(0); i < 8; ++i)
    {
      if(f_pp[i] != f_ptr[i])
      {
        return false;
      }
    }

    // compare index arrays
    for(int i(0); i < 17; ++i)
    {
      if(f_di[i] != f_idx[i])
      {
        return false;
      }
    }

    return true;
  }

  bool test_fg(Graph& fg) const
  {
    CONTEXT("GraphTest::test_fg()");

    // fetch the graph's arrays
    Index* fg_pp = fg.get_domain_ptr();
    Index* fg_di = fg.get_image_idx();

    // check against analytic solution
    //      0  1  2  3  4  5  6
    //   +---------------------
    // 0 |  0  4  1  2  5  3  .
    // 1 |  6  7  8 10  .  9 11
    // 2 | 12 18 13 14 16 15 17
    // 3 | 19 23 20 21 15 22 24
    // 4 | 29  . 26 30 27  . 28
    // 5 | 31 35 32 33  . 34  .
    // 6 |  . 39 36 40 37  . 38

    Index fg_ptr[8] = {0, 6, 12, 19, 26, 31, 36, 41};
    Index fg_idx[41] =
    {
      0, 2, 3, 5, 1, 4,
      0, 1, 2, 5, 3, 6,
      0, 2, 3, 5, 4, 6, 1,
      0, 2, 3, 5, 1, 6, 4,
      2, 4, 6, 0, 3,
      0, 2, 3, 5, 1,
      2, 4, 6, 1, 3
    };

    // check dimensions
    if(fg.get_num_nodes_domain() != 7)
      return false;
    if(fg.get_num_nodes_image() != 7)
      return false;

    // check degree
    if(fg.degree() != 7)
      return false;

    // compare pointer arrays
    for(int i(0); i < 8; ++i)
    {
      if(fg_pp[i] != fg_ptr[i])
      {
        return false;
      }
    }

    // compare index arrays
    for(int i(0); i < 41; ++i)
    {
      if(fg_di[i] != fg_idx[i])
      {
        return false;
      }
    }

    return true;
  }

  bool test_gf(Graph& gf) const
  {
    CONTEXT("GraphTest::test_gf()");

    // fetch the graph's arrays
    Index* gf_pp = gf.get_domain_ptr();
    Index* gf_di = gf.get_image_idx();

    // check against analytic solution
    //      0  1  2  3  4
    //   +---------------
    // 0 |  0  3  1  4  2
    // 1 |  5  6  7  9  8
    // 2 | 10 14 11 13 12
    // 3 | 17 19 15 16 18
    // 4 | 20 24 21 23 22

    Index gf_ptr[6] = {0, 5, 10, 15, 20, 25};
    Index gf_idx[25] =
    {
      0, 2, 4, 1, 3,
      0, 1, 2, 4, 3,
      0, 2, 4, 3, 1,
      2, 3, 0, 4, 1,
      0, 2, 4, 3, 1
    };

    // check dimensions
    if(gf.get_num_nodes_domain() != 5)
      return false;
    if(gf.get_num_nodes_image() != 5)
      return false;

    // check degree
    if(gf.degree() != 5)
      return false;

    // compare pointer arrays
    for(int i(0); i < 6; ++i)
    {
      if(gf_pp[i] != gf_ptr[i])
      {
        return false;
      }
    }

    // compare index arrays
    for(int i(0); i < 25; ++i)
    {
      if(gf_di[i] != gf_idx[i])
      {
        return false;
      }
    }

    return true;
  }

  bool test_constr_perm(Graph& f) const
  {
    CONTEXT("GraphTest::test_constr_perm()");

    // fetch the graph's arrays
    Index* domain_ptr = f.get_domain_ptr();
    Index* image_ptr = f.get_image_idx();

    // check against analytic solution
    //      0  1  2  3  4  5  6
    //   +---------------------
    // 0 |  0  2  3  .  .  .  1
    // 1 |  .  4  .  6  5  .  .
    // 2 |  7  8 10  .  .  9  .
    // 3 |  .  .  . 13  . 12 11
    // 4 | 14  .  .  . 16 15  .

    Index domain_ref[6] = {0, 4, 7, 11, 14, 17};
    Index image_ref[17] =
    {
      0, 6, 1, 2, 1, 4, 3, 0, 1, 5, 2, 6, 5, 3, 0, 5, 4
    };

    // check dimensions
    if(f.get_num_nodes_domain() != 5)
      return false;
    if(f.get_num_nodes_image() != 7)
      return false;

    // compare pointer arrays
    for(int i(0); i < 6; ++i)
    {
      if(domain_ptr[i] != domain_ref[i])
      {
        return false;
      }
    }

    // compare index arrays
    for(int i(0); i < 17; ++i)
    {
      if(image_ptr[i] != image_ref[i])
      {
        return false;
      }
    }
    return true;
  }

  bool test_constr_colour(Graph& f) const
  {
    CONTEXT("GraphTest::test_constr_colour()");

    // fetch the graph's arrays
    Index* domain_ptr = f.get_domain_ptr();
    Index* image_ptr = f.get_image_idx();

    // check against analytic solution
    //      0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
    //   +---------------------------------------------
    // 0 |  .  .  .  .  .  .  .  .  .  .  .  0  1  .  .
    // 1 |  2  .  3  4  .  .  .  .  5  .  .  .  .  6  .
    // 2 |  .  7  .  .  .  .  8  .  .  .  .  .  .  .  9
    // 3 |  .  .  .  . 10  .  . 11  .  .  .  .  .  .  .
    // 4 |  .  .  .  .  . 12  .  .  . 13 14  .  .  .  .

    // corresponding colouring (see run()):
    // nodes:   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
    // colours: 1, 2, 1, 1, 3, 4, 2, 3, 1, 4, 4, 0, 0, 1, 2

    Index domain_ref[6] = {0, 2, 7, 10, 12, 15};
    Index image_ref[15] =
    {
      11, 12,
      0, 2, 3, 8, 13,
      1, 6, 14,
      4, 7,
      5, 9, 10
    };

    // check dimensions
    if(f.get_num_nodes_domain() != 5)
      return false;
    if(f.get_num_nodes_image() != 15)
      return false;

    // compare pointer arrays
    for(int i(0); i < 6; ++i)
    {
      if(domain_ptr[i] != domain_ref[i])
      {
        return false;
      }
    }

    // compare index arrays
    for(int i(0); i < 15; ++i)
    {
      if(image_ptr[i] != image_ref[i])
      {
        return false;
      }
    }
    return true;
  }

  bool test_sort() const
  {
    // create an unsorted graph
    Index g_ptr[8] = {0, 6, 12, 19, 26, 31, 36, 41};
    Index g_idx[41] =
    {
      0, 2, 3, 5, 1, 4,
      0, 1, 2, 5, 3, 6,
      0, 2, 3, 5, 4, 6, 1,
      0, 2, 3, 5, 1, 6, 4,
      2, 4, 6, 0, 3,
      0, 2, 3, 5, 1,
      2, 4, 6, 1, 3
    };

    Graph g(7, 7, 41, g_ptr, nullptr, g_idx);

    // sort the graph
    g.sort_indices();

    // validate the sorted graph's indices
    Index s_idx[41] =
    {
      0, 1, 2, 3, 4, 5,
      0, 1, 2, 3, 5, 6,
      0, 1, 2, 3, 4, 5, 6,
      0, 1, 2, 3, 4, 5, 6,
      0, 2, 3, 4, 6,
      0, 1, 2, 3, 5,
      1, 2, 3, 4, 6
    };
    Index* gsi = g.get_image_idx();
    for(Index i(0); i < 41; ++i)
    {
      if(s_idx[i] != gsi[i])
        return false;
    }

    // okay
    return true;
  }

  virtual void run() const override
  {
    // create a graph G
    //      0  1  2  3  4  5  6
    //   +---------------------
    // 0 |  0  .  1  2  .  3  .
    // 1 |  .  .  4  .  5  .  6
    // 2 |  7  8  9  .  . 10  .
    // 3 |  . 11  . 12  .  . 13
    // 4 | 14  .  . 15 16  .  .

    Index g_ptr[6] = {0, 4, 7, 11, 14, 17};
    Index g_idx[17] =
    {
      0, 2, 3, 5,
      2, 4, 6,
      0, 1, 2, 5,
      1, 3, 6,
      0, 3, 4
    };
    Graph g(5, 7, 17, g_ptr, nullptr, g_idx);

    // transpose the graph G
    Graph f(rt_transpose, g);
    TEST_CHECK(test_f(f));

    // render and test fg
    Graph fg(rt_injectify, f, g);
    TEST_CHECK(test_fg(fg));

    // render and test gf
    Graph gf(rt_injectify, g, f);
    TEST_CHECK(test_gf(gf));

    // test sorting
    TEST_CHECK(test_sort());

    // test the construction with a domain- and an image-permutation

    // define domain permutation
    Index domain_perm_idx[5] =
    {
      2, 1, 0, 3, 4
    };

    Permutation prm_domain(5, Permutation::type_perm, domain_perm_idx);

    // define image permutation
    Index image_perm_idx[7] =
    {
      0, 6, 1, 5, 4, 2, 3
    };

    Permutation prm_image(7, Permutation::type_perm, image_perm_idx);

    // construct graph
    Graph fperm(g, prm_domain, prm_image);

    // test
    TEST_CHECK(test_constr_perm(fperm));

    // test the creation out of a colouring object

    // colouring array
    Index colour[15] =
    {
      1, 2, 1, 1, 3, 4, 2, 3, 1, 4, 4, 0, 0, 1, 2
    };

    // create colouring object
    Colouring col(15, colour);

    // construct graph
    Graph fcol(col.create_partition_graph());

    // validate
    TEST_CHECK(test_constr_colour(fcol));

  }

} graph_test;
