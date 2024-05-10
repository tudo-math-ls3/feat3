// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <stdint.h>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Adjacency;

/**
 * \brief Test class for the Permutation class.
 *
 * \test Tests the Permutation class.
 *
 * \author Peter Zajac
 */
class PermutationTest
  : public UnitTest
{
public:
  PermutationTest() :
    UnitTest("PermutationTest")
  {
  }

  virtual ~PermutationTest()
  {
  }

  template<int n_>
  static void make_id(Index (&v)[n_])
  {
    for(Index i(0); i < n_; ++i)
      v[i] = i;
  }

  template<int n_>
  static bool is_id(const Index (&v)[n_])
  {
    for(Index i(0); i < n_; ++i)
    {
      if(v[i] != i)
        return false;
    }
    return true;
  }

  static bool is_id(const Permutation &p)
  {
      for (Index i(0); i < p.size(); ++i)
      {
          if (p.get_perm_pos()[i] != i)
              return false;
      }
      return true;
  }

  template<int n_>
  static void copy(Index (&v)[n_], const Index (&w)[n_])
  {
    for(Index i(0); i < n_; ++i)
      v[i] = w[i];
  }

  template<int n_>
  static bool is_equal(const Index (&v)[n_], const Index (&w)[n_])
  {
    for(Index i(0); i < n_; ++i)
    {
      if(v[i] != w[i])
        return false;
    }
    return true;
  }

  template<int n_>
  bool test_inv(Permutation& prm_rnd, const Index (&v)[n_]) const
  {
    Index ws[n_], wp[n_];
    copy(ws, v);

    // apply inverse permutation
    prm_rnd.apply(wp, v, true);

    // appyl inverse swapping
    prm_rnd.apply(ws, true);

    // test for id
    return is_id(wp) && is_id(ws);
  }

  template<int n_>
  bool test_fwd(Permutation& prm_inv, const Index (&v)[n_]) const
  {
    Index ws[n_], wp[n_];
    copy(ws, v);

    // apply forward permutation
    prm_inv.apply(wp, v);

    // appyl forward swapping
    prm_inv.apply(ws);

    // test for id
    return is_id(wp) && is_id(ws);
  }

  template<int n_>
  static void make_rand_vec(Index(&v)[n_], Random& rng)
  {
    for (Index i(0); i < n_; ++i)
    {
      v[i] = rng.next();
    }
  }

  virtual void run() const override
  {
    // create an rng
    Random rng;

    static constexpr Index N = 10;
    // create a random permutation
    Permutation prm_rnd0(N, rng);
    Permutation prm_rnd(N, rng);
    prm_rnd = std::move(prm_rnd0);

    // create two identity arrays
    Index vp[N], vs[N];
    make_id(vs);


    // apply permutation operator
    prm_rnd.apply(vp, vs);

    // apply swapping operator
    prm_rnd.apply(vs);

    // make sure they're equal
    TEST_CHECK(is_equal(vp,vs));

    // test inverse permutation
    test_inv(prm_rnd, vs);

    // create inverse random permutation directly
    Permutation prm_inv1(prm_rnd.inverse());
    test_fwd(prm_inv1, vs);

    // create inverse random permutation by inverse permute array
    Permutation prm_inv2(N, Permutation::type_inv_perm, prm_rnd.get_perm_pos());
    test_fwd(prm_inv2, vs);

    // create inverse random permutation by inverse swap array
    Permutation prm_inv3(N, Permutation::type_inv_swap, prm_rnd.get_swap_pos());
    test_fwd(prm_inv3, vs);

    //Test concat
    //create random array
    Index randvec1[N], randvec2[N];
    make_rand_vec(randvec1, rng);
    copy(randvec2, randvec1);

    //create two random permutations
    Permutation prm_rnd1(N, rng);
    Permutation prm_rnd2 (N, rng);

    //first apply prm_rnd2 then prm_rnd1
    prm_rnd2.apply(randvec1);
    prm_rnd1.apply(randvec1);

    //apply the concatenation
    prm_rnd1.concat(prm_rnd2);
    prm_rnd1.apply(randvec2);

    //make sure they're equal
    TEST_CHECK(is_equal(randvec1, randvec2));
  }
} permutation_test;
