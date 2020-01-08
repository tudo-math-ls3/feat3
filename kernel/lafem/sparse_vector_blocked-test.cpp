// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse vector blocked class.
 *
 * \test test description missing
 *
 * \tparam Mem_
 * description missing
 *
 * \tparam DT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseVectorBlockedTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  SparseVectorBlockedTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseVectorBlockedTest")
  {
  }

  virtual ~SparseVectorBlockedTest()
  {
  }

  virtual void run() const override
  {
    SparseVectorBlocked<Mem_, DT_, IT_, 2> zero1;
    SparseVectorBlocked<Mem::Main, DT_, IT_, 2> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    SparseVectorBlocked<Mem_, DT_, IT_, 2> a(10);
    TEST_CHECK_EQUAL(a, a);
    Tiny::Vector<DT_, 2> tv1(41);
    Tiny::Vector<DT_, 2> tv2(42);
    a(3, tv1);
    a(3, tv2);
    a(6, tv1);
    a(3, tv1);
    a(1, tv1);
    a(6, tv2);
    TEST_CHECK_EQUAL(a(3).v[0], tv1.v[0]);
    TEST_CHECK_EQUAL(a(1)[0], tv1[1]);
    TEST_CHECK_EQUAL(a(6)[1], tv2[1]);
    TEST_CHECK_EQUAL(a(2)[0], a.zero_element()[0]);
    TEST_CHECK_EQUAL(a.used_elements(), Index(3));

    SparseVectorBlocked<Mem_, DT_, IT_, 2> b(a.clone());
    TEST_CHECK_EQUAL(a, b);
    a(2, tv2);
    TEST_CHECK_NOT_EQUAL(a, b);

    //increase vector size above alloc_increment
    SparseVectorBlocked<Mem_, DT_, IT_, 2> c(1001);
    for (Index i(1) ; i <= c.size() ; ++i)
    {
      c(c.size() - i, tv1);
    }

    Random::SeedType seed(Random::SeedType(time(nullptr)));
    std::cout << "seed: " << seed << std::endl;
    Random rng(seed);
    Adjacency::Permutation prm_rnd(a.size(), rng);
    SparseVectorBlocked<Mem_, DT_, IT_, 2> ap(a.clone());
    ap.permute(prm_rnd);
    prm_rnd = prm_rnd.inverse();
    ap.permute(prm_rnd);
    TEST_CHECK_EQUAL(ap, a);
    TEST_CHECK_EQUAL(ap.used_elements(), Index(4));
  }
};
SparseVectorBlockedTest<Mem::Main, float, Index> cpu_sparse_vector_blocked_test_float;
SparseVectorBlockedTest<Mem::Main, double, Index> cpu_sparse_vector_blocked_test_double;
#ifdef FEAT_HAVE_CUDA
SparseVectorBlockedTest<Mem::CUDA, float, Index> cuda_sparse_vector_blocked_test_float;
SparseVectorBlockedTest<Mem::CUDA, double, Index> cuda_sparse_vector_blocked_test_double;
#endif
