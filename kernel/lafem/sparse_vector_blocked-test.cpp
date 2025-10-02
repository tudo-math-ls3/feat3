// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>
#include <kernel/util/binary_stream.hpp>

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
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    SparseVectorBlocked<DT_, IT_, 2> zero1;
    SparseVectorBlocked<DT_, IT_, 2> zero2;
    //TEST_CHECK_LESS_THAN(zero1.max_rel_diff(zero2), eps);

    SparseVectorBlocked<DT_, IT_, 2> a(10);
    //TEST_CHECK_LESS_THAN(a.max_rel_diff(a), eps);
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
    TEST_CHECK_LESS_THAN(a.max_rel_diff(b), eps);
    a(3, tv2);
    TEST_CHECK_LESS_THAN(eps, a.max_rel_diff(b));

    //increase vector size above alloc_increment
    SparseVectorBlocked<DT_, IT_, 2> c(1001);
    for (Index i(1) ; i <= c.size() ; ++i)
    {
      c(c.size() - i, tv1);
    }

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";
    Adjacency::Permutation prm_rnd(a.size(), rng);
    SparseVectorBlocked<DT_, IT_, 2> ap(a.clone());
    ap.permute(prm_rnd);
    prm_rnd = prm_rnd.inverse();
    ap.permute(prm_rnd);
    TEST_CHECK_LESS_THAN(ap.max_rel_diff(a), eps);
    TEST_CHECK_EQUAL(ap.used_elements(), Index(3));
  }
};
SparseVectorBlockedTest <float, std::uint32_t> cpu_sparse_vector_blocked_test_float_uint32(PreferredBackend::generic);
SparseVectorBlockedTest <double, std::uint32_t> cpu_sparse_vector_blocked_test_double_uint32(PreferredBackend::generic);
SparseVectorBlockedTest <float, std::uint64_t> cpu_sparse_vector_blocked_test_float_uint64(PreferredBackend::generic);
SparseVectorBlockedTest <double, std::uint64_t> cpu_sparse_vector_blocked_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseVectorBlockedTest <float, std::uint64_t> mkl_cpu_sparse_vector_blocked_test_float_uint64(PreferredBackend::mkl);
SparseVectorBlockedTest <double, std::uint64_t> mkl_cpu_sparse_vector_blocked_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseVectorBlockedTest <__float128, std::uint64_t> cpu_sparse_vector_blocked_test_float128_uint64(PreferredBackend::generic);
SparseVectorBlockedTest <__float128, std::uint32_t> cpu_sparse_vector_blocked_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseVectorBlockedTest <Half, std::uint32_t> cpu_sparse_vector_blocked_test_half_uint32(PreferredBackend::generic);
SparseVectorBlockedTest <Half, std::uint64_t> cpu_sparse_vector_blocked_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseVectorBlockedTest <float, std::uint32_t> cuda_sparse_vector_blocked_test_float_uint32(PreferredBackend::cuda);
SparseVectorBlockedTest <double, std::uint32_t> cuda_sparse_vector_blocked_test_double_uint32(PreferredBackend::cuda);
SparseVectorBlockedTest <float, std::uint64_t> cuda_sparse_vector_blocked_test_float_uint64(PreferredBackend::cuda);
SparseVectorBlockedTest <double, std::uint64_t> cuda_sparse_vector_blocked_test_double_uint64(PreferredBackend::cuda);
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
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    SparseVectorBlocked<DT_, IT_, 2> a(10);
    //TEST_CHECK_LESS_THAN(a.max_rel_diff(a), eps);
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
    TEST_CHECK_LESS_THAN(bin.max_rel_diff(a), eps);

    auto op = a.serialize(LAFEM::SerialConfig(false, false));
    SparseVectorBlocked<DT_, IT_, 2> o(op);
    TEST_CHECK_LESS_THAN(a.max_rel_diff(o), eps);
#ifdef FEAT_HAVE_ZLIB
    auto zl = a.serialize(LAFEM::SerialConfig(true, false));
    SparseVectorBlocked<DT_, IT_, 2> zlib(zl);
    TEST_CHECK_LESS_THAN(zlib.max_rel_diff(a), eps);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = a.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseVectorBlocked<DT_, IT_, 2> zfp(zf);
    for (Index i(0) ; i < a.size() ; ++i)
    {
      for(int j(0) ; j < a(i).n ; ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(zfp(i)[j], a(i)[j], Math::pow(Math::eps<DT_>(), DT_(0.7)));
      }
    }
#endif
  }
};
SparseVectorBlockedSerializeTest <float, std::uint32_t> cpu_sparse_vector_blocked_serialize_test_float_uint32(PreferredBackend::generic);
SparseVectorBlockedSerializeTest <double, std::uint32_t> cpu_sparse_vector_blocked_serialize_test_double_uint32(PreferredBackend::generic);
SparseVectorBlockedSerializeTest <float, std::uint64_t> cpu_sparse_vector_blocked_serialize_test_float_uint64(PreferredBackend::generic);
SparseVectorBlockedSerializeTest <double, std::uint64_t> cpu_sparse_vector_blocked_serialize_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseVectorBlockedSerializeTest <float, std::uint64_t> mkl_cpu_sparse_vector_blocked_serialize_test_float_uint64(PreferredBackend::mkl);
SparseVectorBlockedSerializeTest <double, std::uint64_t> mkl_cpu_sparse_vector_blocked_serialize_test_double_uint64(PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SparseVectorBlockedSerializeTest <__float128, std::uint64_t> cpu_sparse_vector_blocked_serialize_test_float128_uint64(PreferredBackend::generic);
//SparseVectorBlockedSerializeTest <__float128, std::uint32_t> cpu_sparse_vector_blocked_serialize_test_float128_uint32(PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SparseVectorBlockedSerializeTest <Half, std::uint32_t> cpu_sparse_vector_blocked_serialize_test_half_uint32(PreferredBackend::generic);
SparseVectorBlockedSerializeTest <Half, std::uint64_t> cpu_sparse_vector_blocked_serialize_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseVectorBlockedSerializeTest <float, std::uint32_t> cuda_sparse_vector_blocked_serialize_test_float_uint32(PreferredBackend::cuda);
SparseVectorBlockedSerializeTest <double, std::uint32_t> cuda_sparse_vector_blocked_serialize_test_double_uint32(PreferredBackend::cuda);
SparseVectorBlockedSerializeTest <float, std::uint64_t> cuda_sparse_vector_blocked_serialize_test_float_uint64(PreferredBackend::cuda);
SparseVectorBlockedSerializeTest <double, std::uint64_t> cuda_sparse_vector_blocked_serialize_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseVectorBlockedMaxRelDiffTest
  : public UnitTest
{
public:
  SparseVectorBlockedMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("SparseVectorBlockedMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseVectorBlockedMaxRelDiffTest()
  {
  }

  virtual void run() const override
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    const DT_ delta = DT_(123.5);
    const DT_ initial_value = DT_(10.0);
    static constexpr int block_size = 2;

    const Index size = 100;
    const Index diff_index = 42;

    // reference vector
    SparseVectorBlocked<DT_, IT_, block_size> b(size);
    Tiny::Vector<DT_, block_size> initial_block(initial_value);

    // b(42) = 10
    b(diff_index, initial_block);

    // copy b into a
    SparseVectorBlocked<DT_, IT_, block_size> a = b.clone();

    // a(diff_index) = initial_value + delta
    Tiny::Vector<DT_, block_size> delta_block(delta);
    a(diff_index, a(diff_index)+delta_block);

    // reference value
    const DT_ ref = delta / (DT_(2) * initial_value + delta);

    // test ||a-b||_infty
    const DT_ diff_1 = a.max_rel_diff(b);
    TEST_CHECK_RELATIVE(diff_1, ref, eps);

    // test ||b-a||_infty
    const DT_ diff_2 = b.max_rel_diff(a);
    TEST_CHECK_RELATIVE(diff_2, ref, eps);
  }
};
SparseVectorBlockedMaxRelDiffTest <float, std::uint32_t> cpu_sparse_vector_blocked_max_rel_diff_test_float_uint32(PreferredBackend::generic);
SparseVectorBlockedMaxRelDiffTest <double, std::uint32_t> cpu_sparse_vector_blocked_max_rel_diff_test_double_uint32(PreferredBackend::generic);
SparseVectorBlockedMaxRelDiffTest <float, std::uint64_t> cpu_sparse_vector_blocked_max_rel_diff_test_float_uint64(PreferredBackend::generic);
SparseVectorBlockedMaxRelDiffTest <double, std::uint64_t> cpu_sparse_vector_blocked_max_rel_diff_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseVectorBlockedMaxRelDiffTest <float, std::uint64_t> mkl_cpu_sparse_vector_blocked_max_rel_diff_test_float_uint64(PreferredBackend::mkl);
SparseVectorBlockedMaxRelDiffTest <double, std::uint64_t> mkl_cpu_sparse_vector_blocked_max_rel_diff_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseVectorBlockedMaxRelDiffTest <__float128, std::uint64_t> cpu_sparse_vector_blocked_max_rel_diff_test_float128_uint64(PreferredBackend::generic);
SparseVectorBlockedMaxRelDiffTest <__float128, std::uint32_t> cpu_sparse_vector_blocked_max_rel_diff_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseVectorBlockedMaxRelDiffTest <Half, std::uint32_t> cpu_sparse_vector_blocked_max_rel_diff_test_half_uint32(PreferredBackend::generic);
SparseVectorBlockedMaxRelDiffTest <Half, std::uint64_t> cpu_sparse_vector_blocked_max_rel_diff_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SparseVectorBlockedMaxRelDiffTest <Half, std::uint32_t> cuda_sparse_vector_blocked_max_rel_diff_test_half_uint32(PreferredBackend::cuda);
SparseVectorBlockedMaxRelDiffTest <Half, std::uint64_t> cuda_sparse_vector_blocked_max_rel_diff_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SparseVectorBlockedMaxRelDiffTest <float, std::uint32_t> cuda_sparse_vector_blocked_max_rel_diff_test_float_uint32(PreferredBackend::cuda);
SparseVectorBlockedMaxRelDiffTest <double, std::uint32_t> cuda_sparse_vector_blocked_max_rel_diff_test_double_uint32(PreferredBackend::cuda);
SparseVectorBlockedMaxRelDiffTest <float, std::uint64_t> cuda_sparse_vector_blocked_max_rel_diff_test_float_uint64(PreferredBackend::cuda);
SparseVectorBlockedMaxRelDiffTest <double, std::uint64_t> cuda_sparse_vector_blocked_max_rel_diff_test_double_uint64(PreferredBackend::cuda);
#endif
