#include <kernel/base_header.hpp>
#include <kernel/util/graph.hpp>
#include <test_system/test_system.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

typedef CompositeAdjactor<Graph,Graph> CAGG;

/**
 * \brief Test class for the Graph class.
 *
 * \test Tests the Graph class.
 *
 * \author Peter Zajac
 */
class GraphTest
  : public TaggedTest<Nil, Nil>
{
public:
  GraphTest() :
    TaggedTest<Nil, Nil>("graph_test")
  {
  }

  bool test_f(Graph* f) const
  {
    CONTEXT("GraphTest::test_f()");

    // fetch the graph's arrays
    Index* f_pp = f->get_domain_ptr();
    Index* f_di = f->get_image_idx();

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

    // check degree
    if(f->degree() != 3)
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

  bool test_fg(Graph* fg) const
  {
    CONTEXT("GraphTest::test_fg()");

    // fetch the graph's arrays
    Index* fg_pp = fg->get_domain_ptr();
    Index* fg_di = fg->get_image_idx();

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

    // check degree
    if(fg->degree() != 7)
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

  bool test_gf(Graph* gf) const
  {
    CONTEXT("GraphTest::test_gf()");

    // fetch the graph's arrays
    Index* gf_pp = gf->get_domain_ptr();
    Index* gf_di = gf->get_image_idx();

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

    // check degree
    if(gf->degree() != 5)
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

  virtual void run() const
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
    Graph* g = new Graph(5, 7, 17, g_ptr, nullptr, g_idx);

    // transpose the graph G
    Graph* f = Graph::render(*g, false, true);
    TEST_CHECK(test_f(f));

    // render and test fg
    Graph* fg = Graph::render_composite(*f, *g, true, false);
    TEST_CHECK(test_fg(fg));

    // render and test gf
    Graph* gf = Graph::render_composite(*g, *f, true, false);
    TEST_CHECK(test_gf(gf));

    // clean up the mess
    delete gf;
    delete fg;
    delete f;
    delete g;
  }

} graph_test;
