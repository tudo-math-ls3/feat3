// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/colouring.hpp>
#include <kernel/adjacency/permutation.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Adjacency;

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

  virtual ~GraphTest()
  {
  }

  static void shuffle_indices(Graph& graph, Random& rng)
  {
    const Index n = graph.get_num_nodes_domain();
    const Index* ptr = graph.get_domain_ptr();
    Index* idx = graph.get_image_idx();

    for(Index i(0); i < n; ++i)
    {
      if(ptr[i+1] <= ptr[i]+1u)
        continue;

      for(Index j(ptr[i]); j+1u < ptr[i+1]; ++j)
      {
        Index k = rng(j, ptr[i+1u]-1u);
        if(k > j)
        {
          Index t = idx[j];
          idx[j] = idx[k];
          idx[k] = t;
        }
      }
    }
  }

  void test_f(Graph& f) const
  {
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

    const Index f_ptr[8] = {0, 3, 5, 8, 11, 13, 15, 17};
    const Index f_idx[17] =
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
    TEST_CHECK_EQUAL(f.get_num_nodes_domain(), 7);
    TEST_CHECK_EQUAL(f.get_num_nodes_image(), 5);

    // check degree
    TEST_CHECK_EQUAL(f.degree(), 3);

    // compare pointer arrays
    for(int i(0); i < 8; ++i)
    {
      TEST_CHECK_EQUAL(f_pp[i], f_ptr[i]);
    }

    // compare index arrays
    for(int i(0); i < 17; ++i)
    {
      TEST_CHECK_EQUAL(f_di[i], f_idx[i]);
    }
  }

  void test_fg(Graph& fg) const
  {
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

    const Index fg_ptr[8] = {0, 6, 12, 19, 26, 31, 36, 41};
    const Index fg_idx[41] =
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
    TEST_CHECK_EQUAL(fg.get_num_nodes_domain(), 7);
    TEST_CHECK_EQUAL(fg.get_num_nodes_image(), 7);

    // check degree
    TEST_CHECK_EQUAL(fg.degree(), 7);

    // compare pointer arrays
    for(int i(0); i < 8; ++i)
    {
      TEST_CHECK_EQUAL(fg_pp[i], fg_ptr[i]);
    }

    // compare index arrays
    for(int i(0); i < 41; ++i)
    {
      TEST_CHECK_EQUAL(fg_di[i], fg_idx[i]);
    }
  }

  void test_gf(Graph& gf) const
  {
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

    const Index gf_ptr[6] = {0, 5, 10, 15, 20, 25};
    const Index gf_idx[25] =
    {
      0, 2, 4, 1, 3,
      0, 1, 2, 4, 3,
      0, 2, 4, 3, 1,
      2, 3, 0, 4, 1,
      0, 2, 4, 3, 1
    };

    // check dimensions
    TEST_CHECK_EQUAL(gf.get_num_nodes_domain(), 5);
    TEST_CHECK_EQUAL(gf.get_num_nodes_image(), 5);

    // check degree
    TEST_CHECK_EQUAL(gf.degree(), 5);

    // compare pointer arrays
    for(int i(0); i < 6; ++i)
    {
      TEST_CHECK_EQUAL(gf_pp[i], gf_ptr[i]);
    }

    // compare index arrays
    for(int i(0); i < 25; ++i)
    {
      TEST_CHECK_EQUAL(gf_di[i], gf_idx[i]);
    }
  }

  void test_basic() const
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
    Graph g(5, 7, 17, g_ptr, g_idx);

    // transpose the graph G
    Graph f(RenderType::transpose, g);
    test_f(f);

    // render and test fg
    Graph fg(RenderType::injectify, f, g);
    test_fg(fg);

    // render and test gf
    Graph gf(RenderType::injectify, g, f);
    test_gf(gf);
  }

  void test_constr_perm() const
  {
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
    Graph g(5, 7, 17, g_ptr, g_idx);

    // construct graph
    Graph f(g, prm_domain, prm_image);

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

    const Index domain_ref[6] = {0, 4, 7, 11, 14, 17};
    const Index image_ref[17] =
    {
      0, 6, 1, 2, 1, 4, 3, 0, 1, 5, 2, 6, 5, 3, 0, 5, 4
    };

    // check dimensions
    TEST_CHECK_EQUAL(f.get_num_nodes_domain(), 5);
    TEST_CHECK_EQUAL(f.get_num_nodes_image(), 7);

    // compare pointer arrays
    for(int i(0); i < 6; ++i)
    {
      TEST_CHECK_EQUAL(domain_ptr[i], domain_ref[i]);
    }

    // compare index arrays
    for(int i(0); i < 17; ++i)
    {
      TEST_CHECK_EQUAL(image_ptr[i], image_ref[i]);
    }
  }

  void test_constr_colour() const
  {
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

    // fetch the graph's arrays
    Index* domain_ptr = fcol.get_domain_ptr();
    Index* image_ptr = fcol.get_image_idx();

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

    const Index domain_ref[6] = {0, 2, 7, 10, 12, 15};
    const Index image_ref[15] =
    {
      11, 12,
      0, 2, 3, 8, 13,
      1, 6, 14,
      4, 7,
      5, 9, 10
    };

    // check dimensions
    TEST_CHECK_EQUAL(fcol.get_num_nodes_domain(), 5);
    TEST_CHECK_EQUAL(fcol.get_num_nodes_image(), 15);

    // compare pointer arrays
    for(int i(0); i < 6; ++i)
    {
      TEST_CHECK_EQUAL(domain_ptr[i], domain_ref[i]);
    }

    // compare index arrays
    for(int i(0); i < 15; ++i)
    {
      TEST_CHECK_EQUAL(image_ptr[i], image_ref[i]);
    }
  }

  void test_sort() const
  {
    // create an unsorted graph
    const Index g_ptr[8] = {0, 6, 12, 19, 26, 31, 36, 41};
    const Index g_idx[41] =
    {
      0, 2, 3, 5, 1, 4,
      0, 1, 2, 5, 3, 6,
      0, 2, 3, 5, 4, 6, 1,
      0, 2, 3, 5, 1, 6, 4,
      2, 4, 6, 0, 3,
      0, 2, 3, 5, 1,
      2, 4, 6, 1, 3
    };

    Graph g(7, 7, 41, g_ptr, g_idx);

    // sort the graph
    g.sort_indices();

    // validate the sorted graph's indices
    const Index s_idx[41] =
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
      TEST_CHECK_EQUAL(s_idx[i], gsi[i]);
    }
  }

  void check_graph(const Graph& g, const Index num_dom, const Index num_img, const Index num_idx,
    const Index* dom_ptr, const Index* img_idx) const
  {
    TEST_CHECK_EQUAL(g.get_num_nodes_domain(), num_dom);
    TEST_CHECK_EQUAL(g.get_num_nodes_image(), num_img);
    TEST_CHECK_EQUAL(g.get_num_indices(), num_idx);

    // get arrays
    const Index* g_dom_ptr = g.get_domain_ptr();
    const Index* g_img_idx = g.get_image_idx();

    // compare pointer array
    for(Index i(0); i <= num_dom; ++i)
    {
      TEST_CHECK_EQUAL(g_dom_ptr[i], dom_ptr[i]);
    }

    // compare index array
    for(Index i(0); i < num_idx; ++i)
    {
      TEST_CHECK_EQUAL(g_img_idx[i], img_idx[i]);
    }
  }

  void check_graph(const Graph& g, const Graph& f) const
  {
    TEST_CHECK_EQUAL(g.get_num_nodes_domain(), f.get_num_nodes_domain());
    TEST_CHECK_EQUAL(g.get_num_nodes_image(), f.get_num_nodes_image());
    TEST_CHECK_EQUAL(g.get_num_indices(), f.get_num_indices());

    const Index num_dom = g.get_num_nodes_domain();
    const Index num_idx = g.get_num_indices();

    // get arrays
    const Index* g_dom_ptr = g.get_domain_ptr();
    const Index* g_img_idx = g.get_image_idx();
    const Index* f_dom_ptr = f.get_domain_ptr();
    const Index* f_img_idx = f.get_image_idx();

    // compare pointer array
    for(Index i(0); i <= num_dom; ++i)
    {
      TEST_CHECK_EQUAL(g_dom_ptr[i], f_dom_ptr[i]);
    }

    // compare index array
    for(Index i(0); i < num_idx; ++i)
    {
      TEST_CHECK_EQUAL(g_img_idx[i], f_img_idx[i]);
    }
  }

  void test_mesh() const
  {
    //   6-----9-----7----11-----8
    //   |           |'\.        |
    //   |           |  '\.   5  |
    //  10     3    13    '8.   12
    //   |           |  4   '\.  |
    //   |           |        '\.|
    //   3-----7-----4-----6-----5
    //   |        ./'|           |
    //   |  0   ./'  |           |
    //   4    .1'    5     2     2
    //   |  ./'   1  |           |
    //   |./'        |           |
    //   0-----0-----1-----3-----2

    // verts: 9
    // edges: 14
    // cells: 6

    // Here come all basic graphs, although not all of them are used.

    // vertices-at-edge
    const Index ve_ptr[15] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28};
    const Index ve_idx[28] = //sorted
    {
      0, 1, // 0
      0, 4,
      2, 5,
      1, 2,
      0, 3,
      1, 4, // 5
      4, 5,
      3, 4,
      5, 7,
      6, 7,
      3, 6, // 10
      7, 8,
      5, 8,
      4, 7
    };

    // vertices-at-cell
    const Index vc_ptr[7] = {0, 3, 6, 10, 14, 17, 20};
    const Index vc_idx[20] = // sorted
    {
      0, 3, 4,
      0, 1, 4,
      1, 2, 4, 5,
      3, 4, 6, 7,
      4, 5, 7,
      5, 7, 8
    };

    // edges-at-vert
    const Index ev_ptr[10] = {0, 3, 6, 8, 11, 15, 19, 21, 24, 26};
    const Index ev_idx[26] = // sorted
    {
      0, 1, 4,
      0, 3, 5,
      2, 3,
      4, 7, 10,
      5, 6, 7, 13,
      2, 6, 8, 12,
      9, 10,
      9, 11, 13,
      11, 12
    };

    // edges-at-cell
    const Index ec_ptr[7] = {0, 3, 6, 10, 14, 17, 20};
    const Index ec_idx[20] = // sorted
    {
      1, 4, 7,
      0, 1, 5,
      2, 3, 5, 6,
      7, 9, 10, 13,
      6, 8, 13,
      8, 11, 12
    };

    // cells-at-vert
    const Index cv_ptr[10] = {0, 2, 4, 5, 7, 12, 15, 16, 19, 20};
    const Index cv_idx[20] = // sorted
    {
      0, 1,
      1, 2,
      2,
      0, 3,
      0, 1, 2, 3, 4,
      2, 4, 5,
      3,
      3, 4, 5,
      5
    };

    // cells-at-edge
    const Index ce_ptr[15] = {0, 1, 3, 4, 5, 6, 8, 10, 12, 14, 15, 16, 17, 18, 20};
    const Index ce_idx[20] = // sorted
    {
      1,    // 0
      0, 1,
      2,
      2,
      0,
      1, 2, // 5
      2, 4,
      0, 3,
      4, 5,
      3,
      3,    // 10
      5,
      5,
      3, 4
    };

    // vertices-over-cell-at-edge = verts-at-cell * cells-at-edge
    const Index vce_ptr[15] = {0, 3, 9, 13, 17, 20, 27, 34, 41, 47, 51, 55, 58, 61, 68};
    const Index vce_idx[68] = // sorted
    {
      0, 1, 4,             // 0
      0, 0, 1, 3, 4, 4,    // 1
      1, 2, 4, 5,          // 2
      1, 2, 4, 5,          // 3
      0, 3, 4,             // 4
      0, 1, 1, 2, 4, 4, 5, // 5
      1, 2, 4, 4, 5, 5, 7, // 6
      0, 3, 3, 4, 4, 6, 7, // 7
      4, 5, 5, 7, 7, 8,    // 8
      3, 4, 6, 7,          // 9
      3, 4, 6, 7,          // 10
      5, 7, 8,             // 11
      5, 7, 8,             // 12
      3, 4, 4, 5, 6, 7, 7  // 13
    };

    // vertices-over-cell-at-edge = verts-at-cell * cells-at-edge
    const Index vce_i_ptr[15] = {0, 3, 7, 11, 15, 18, 23, 28, 33, 37, 41, 45, 48, 51, 56};
    const Index vce_i_idx[56] = // sorted + injectified
    {
      0, 1, 4,       // 0
      0, 1, 3, 4,    // 1
      1, 2, 4, 5,    // 2
      1, 2, 4, 5,    // 3
      0, 3, 4,       // 4
      0, 1, 2, 4, 5, // 5
      1, 2, 4, 5, 7, // 6
      0, 3, 4, 6, 7, // 7
      4, 5, 7, 8,    // 8
      3, 4, 6, 7,    // 9
      3, 4, 6, 7,    // 10
      5, 7, 8,       // 11
      5, 7, 8,       // 12
      3, 4, 5, 6, 7  // 13
    };

    // edges-over-cell-at-vertex = edges-at-cell * cells-at-vert
    const Index ecv_ptr[10] = {0, 6, 13, 17, 24, 41, 51, 55, 65, 68};
    const Index ecv_idx[68] = // sorted
    {
      0, 1, 1, 4, 5, 7,                  // 0
      0, 1, 2, 3, 5, 5, 6,               // 1
      2, 3, 5, 6,                        // 2
      1, 4, 7, 7, 9, 10, 13,             // 3
      0, 1, 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 9, 10, 13, 13, // 4
      2, 3, 5, 6, 6, 8, 8, 11, 12, 13,   // 5
      7, 9, 10, 13,                      // 6
      6, 7, 8, 8, 9, 10, 11, 12, 13, 13, // 7
      8, 11, 12                          // 8
    };

    // edges-over-cell-at-vertex = edges-at-cell * cells-at-vert
    const Index ecv_i_ptr[10] = {0, 5, 11, 15, 21, 33, 41, 45, 53, 56};
    const Index ecv_i_idx[56] = // sorted + injectified
    {
      0, 1, 4, 5, 7,                        // 0
      0, 1, 2, 3, 5, 6,                     // 1
      2, 3, 5, 6,                           // 2
      1, 4, 7, 9, 10, 13,                   // 3
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, // 4
      2, 3, 5, 6, 8, 11, 12, 13,            // 5
      7, 9, 10, 13,                         // 6
      6, 7, 8, 9, 10, 11, 12, 13,           // 7
      8, 11, 12                             // 8
    };

    // create graphs and shuffle them
    Random rng;
    Graph g_ve(14, 9, 28, ve_ptr, ve_idx);    shuffle_indices(g_ve, rng);
    Graph g_vc(6, 9, 20, vc_ptr, vc_idx);     shuffle_indices(g_vc, rng);
    Graph g_ev(9, 14, 26, ev_ptr, ev_idx);    shuffle_indices(g_ev, rng);
    Graph g_ec(6, 14, 20, ec_ptr, ec_idx);    shuffle_indices(g_ec, rng);
    Graph g_cv(9, 6, 20, cv_ptr, cv_idx);     shuffle_indices(g_cv, rng);
    Graph g_ce(14, 6, 20, ce_ptr, ce_idx);    shuffle_indices(g_ce, rng);
    Graph g_vce(14, 9, 68, vce_ptr, vce_idx); shuffle_indices(g_vce, rng);
    Graph g_ecv(9, 14, 68, ecv_ptr, ecv_idx); shuffle_indices(g_ecv, rng);

    Graph graph;

    // test as-is
    graph = Graph(RenderType::as_is, g_cv);
    check_graph(graph, g_cv);

    // test as_is_sorted
    graph = Graph(RenderType::as_is_sorted, g_cv);
    check_graph(graph, 9, 6, 20, cv_ptr, cv_idx);

    // test injectify_sorted
    graph = Graph(RenderType::injectify_sorted, g_vce);
    check_graph(graph, 14, 9, 56, vce_i_ptr, vce_i_idx);

    // test transpose (auto-sorted)
    graph = Graph(RenderType::transpose, g_ce);
    check_graph(graph, 6, 14, 20, ec_ptr, ec_idx);

    // test injectify_transpose (auto-sorted)
    graph = Graph(RenderType::injectify_transpose, g_vce);
    check_graph(graph, 9, 14, 56, ecv_i_ptr, ecv_i_idx);

    // test composite as_is_sorted
    graph = Graph(RenderType::as_is_sorted, g_ce, g_vc);
    check_graph(graph, 14, 9, 68, vce_ptr, vce_idx);

    // test composite injectify_sorted
    graph = Graph(RenderType::injectify_sorted, g_ce, g_vc);
    check_graph(graph, 14, 9, 56, vce_i_ptr, vce_i_idx);

    // test composite transpose (auto-sorted)
    graph = Graph(RenderType::transpose, g_ce, g_vc);
    check_graph(graph, 9, 14, 68, ecv_ptr, ecv_idx);

    // test composite injectify_transpose (auto-sorted)
    graph = Graph(RenderType::injectify_transpose, g_ce, g_vc);
    check_graph(graph, 9, 14, 56, ecv_i_ptr, ecv_i_idx);
  }

  virtual void run() const override
  {
    // perform basic
    test_basic();

    // test sorting
    test_sort();

    // test mesh
    test_mesh();

    // test permutation constructor
    test_constr_perm();

    // test colour constructor
    test_constr_colour();
  }

} graph_test;
