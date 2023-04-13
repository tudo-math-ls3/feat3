// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse vector blocked class.
 *
 * \test test description missing
 *
 * \tparam DT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseVectorBlockedTest
  : public UnitTest
{
public:
  SparseVectorBlockedTest(PreferredBackend backend)
    : UnitTest("SparseVectorBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseVectorBlockedTest()
  {
  }

  virtual void run() const override
  {
    SparseVectorBlocked<DT_, IT_, 2> zero1;
    SparseVectorBlocked<DT_, IT_, 2> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    SparseVectorBlocked<DT_, IT_, 2> a(10);
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
    TEST_CHECK_EQUAL(a.used_elements(), Index(3));

    SparseVectorBlocked<DT_, IT_, 2> b(a.clone());
    TEST_CHECK_EQUAL(a, b);
    a(2, tv2);
    TEST_CHECK_NOT_EQUAL(a, b);

    //increase vector size above alloc_increment
    SparseVectorBlocked<DT_, IT_, 2> c(1001);
    for (Index i(1) ; i <= c.size() ; ++i)
    {
      c(c.size() - i, tv1);
    }

    Random::SeedType seed(Random::SeedType(time(nullptr)));
    std::cout << "seed: " << seed << std::endl;
    Random rng(seed);
    Adjacency::Permutation prm_rnd(a.size(), rng);
    SparseVectorBlocked<DT_, IT_, 2> ap(a.clone());
    ap.permute(prm_rnd);
    prm_rnd = prm_rnd.inverse();
    ap.permute(prm_rnd);
    TEST_CHECK_EQUAL(ap, a);
    TEST_CHECK_EQUAL(ap.used_elements(), Index(4));
  }
};
SparseVectorBlockedTest<float, unsigned int> cpu_sparse_vector_blocked_test_float_uint(PreferredBackend::generic);
SparseVectorBlockedTest<double, unsigned int> cpu_sparse_vector_blocked_test_double_uint(PreferredBackend::generic);
SparseVectorBlockedTest<float, unsigned long> cpu_sparse_vector_blocked_test_float_ulong(PreferredBackend::generic);
SparseVectorBlockedTest<double, unsigned long> cpu_sparse_vector_blocked_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseVectorBlockedTest<float, unsigned long> mkl_cpu_sparse_vector_blocked_test_float_ulong(PreferredBackend::mkl);
SparseVectorBlockedTest<double, unsigned long> mkl_cpu_sparse_vector_blocked_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseVectorBlockedTest<__float128, unsigned long> cpu_sparse_vector_blocked_test_float128_ulong(PreferredBackend::generic);
SparseVectorBlockedTest<__float128, unsigned int> cpu_sparse_vector_blocked_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseVectorBlockedTest<Half, unsigned int> cpu_sparse_vector_blocked_test_half_uint(PreferredBackend::generic);
SparseVectorBlockedTest<Half, unsigned long> cpu_sparse_vector_blocked_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseVectorBlockedTest<float, unsigned int> cuda_sparse_vector_blocked_test_float_uint(PreferredBackend::cuda);
SparseVectorBlockedTest<double, unsigned int> cuda_sparse_vector_blocked_test_double_uint(PreferredBackend::cuda);
SparseVectorBlockedTest<float, unsigned long> cuda_sparse_vector_blocked_test_float_ulong(PreferredBackend::cuda);
SparseVectorBlockedTest<double, unsigned long> cuda_sparse_vector_blocked_test_double_ulong(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseVectorBlockedSerializeTest
  : public UnitTest
{
public:
  SparseVectorBlockedSerializeTest(PreferredBackend backend)
    : UnitTest("SparseVectorBlockedSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseVectorBlockedSerializeTest()
  {
  }

  virtual void run() const override
  {
    SparseVectorBlocked<DT_, IT_, 2> a(10);
    TEST_CHECK_EQUAL(a, a);
    Tiny::Vector<DT_, 2> tv1(41);
    Tiny::Vector<DT_, 2> tv2(42);
    a(3, tv1);
    a(3, tv2);
    a(6, tv1);
    a(3, tv1);
    a(1, tv1);
    a(6, tv2);

    BinaryStream bs;
    a.write_out(FileMode::fm_svb, bs);
    bs.seekg(0);
    SparseVectorBlocked<DT_, IT_, 2> bin(FileMode::fm_svb, bs);
    TEST_CHECK_EQUAL(bin, a);

    auto op = a.serialize(LAFEM::SerialConfig(false, false));
    SparseVectorBlocked<DT_, IT_, 2> o(op);
    TEST_CHECK_EQUAL(a, o);
#ifdef FEAT_HAVE_ZLIB
    auto zl = a.serialize(LAFEM::SerialConfig(true, false));
    SparseVectorBlocked<DT_, IT_, 2> zlib(zl);
    TEST_CHECK_EQUAL(zlib, a);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = a.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseVectorBlocked<DT_, IT_, 2> zfp(zf);
    for (Index i(0) ; i < a.size() ; ++i)
    {
      for(int j(0) ; j < a(i).n ; ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(zfp(i)[j], a(i)[j], DT_(1e-4));
      }
    }
#endif
  }
};
SparseVectorBlockedSerializeTest<float, unsigned int> cpu_sparse_vector_blocked_serialize_test_float_uint(PreferredBackend::generic);
SparseVectorBlockedSerializeTest<double, unsigned int> cpu_sparse_vector_blocked_serialize_test_double_uint(PreferredBackend::generic);
SparseVectorBlockedSerializeTest<float, unsigned long> cpu_sparse_vector_blocked_serialize_test_float_ulong(PreferredBackend::generic);
SparseVectorBlockedSerializeTest<double, unsigned long> cpu_sparse_vector_blocked_serialize_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseVectorBlockedSerializeTest<float, unsigned long> mkl_cpu_sparse_vector_blocked_serialize_test_float_ulong(PreferredBackend::mkl);
SparseVectorBlockedSerializeTest<double, unsigned long> mkl_cpu_sparse_vector_blocked_serialize_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseVectorBlockedSerializeTest<__float128, unsigned long> cpu_sparse_vector_blocked_serialize_test_float128_ulong(PreferredBackend::generic);
SparseVectorBlockedSerializeTest<__float128, unsigned int> cpu_sparse_vector_blocked_serialize_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseVectorBlockedSerializeTest<Half, unsigned int> cpu_sparse_vector_blocked_serialize_test_half_uint(PreferredBackend::generic);
SparseVectorBlockedSerializeTest<Half, unsigned long> cpu_sparse_vector_blocked_serialize_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseVectorBlockedSerializeTest<float, unsigned int> cuda_sparse_vector_blocked_serialize_test_float_uint(PreferredBackend::cuda);
SparseVectorBlockedSerializeTest<double, unsigned int> cuda_sparse_vector_blocked_serialize_test_double_uint(PreferredBackend::cuda);
SparseVectorBlockedSerializeTest<float, unsigned long> cuda_sparse_vector_blocked_serialize_test_float_ulong(PreferredBackend::cuda);
SparseVectorBlockedSerializeTest<double, unsigned long> cuda_sparse_vector_blocked_serialize_test_double_ulong(PreferredBackend::cuda);
#endif
