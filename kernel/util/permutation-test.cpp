#include <test_system/test_system.hpp>
#include <kernel/util/permutation.hpp>
#include <iostream>
#include <ctime>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for the Permutation class.
 *
 * \test Tests the Permutation class.
 *
 * \author Peter Zajac
 */
class PermutationTest
  : public TaggedTest<Archs::None, Nil>
{
public:
  PermutationTest() :
    TaggedTest<Archs::None, Nil>("PermutationTest")
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
    prm_rnd(wp, v, true);

    // appyl inverse swapping
    prm_rnd(ws, true);

    // test for id
    return is_id(wp) && is_id(ws);
  }

  template<int n_>
  bool test_fwd(Permutation& prm_inv, const Index (&v)[n_]) const
  {
    Index ws[n_], wp[n_];
    copy(ws, v);

    // apply forward permutation
    prm_inv(wp, v);

    // appyl forward swapping
    prm_inv(ws);

    // test for id
    return is_id(wp) && is_id(ws);
  }

  virtual void run() const
  {
    // get a seed for us and print it
    uint32_t seed((uint32_t)time(NULL));
    std::cout << "seed = " << seed << std::endl;

    // create an rng
    Random rng(seed);

#define N 10
    // create a random permutation
    Permutation prm_rnd(N, rng);

    // create two identity arrays
    Index vp[N], vs[N];
    make_id(vs);

    // apply permutation operator
    prm_rnd(vp, vs);

    // apply swapping operator
    prm_rnd(vs);

    // make sure they're equal
    TEST_CHECK(is_equal(vp,vs));

    // test inverse permutation
    test_inv(prm_rnd, vs);

    // create inverse random permutation directory
    Permutation prm_inv1(prm_rnd, true);
    test_fwd(prm_inv1, vs);

    // create inverse random permutation by inverse permute array
    Permutation prm_inv2(N, Permutation::type_inv_perm, prm_rnd.get_perm_pos());
    test_fwd(prm_inv2, vs);

    // create inverse random permutation by inverse swap array
    Permutation prm_inv3(N, Permutation::type_inv_swap, prm_rnd.get_swap_pos());
    test_fwd(prm_inv3, vs);
#undef N
  }
} permutation_test;
