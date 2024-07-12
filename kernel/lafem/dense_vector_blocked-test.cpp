// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/binary_stream.hpp>

#include <list>
#include <sstream>
#include <cstdio>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the dense vector blocked class.
 *
 * \test test description missing
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class DenseVectorBlockedTest
  : public UnitTest
{
public:
  DenseVectorBlockedTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedTest()
  {
  }

  virtual void run() const override
  {
    DenseVectorBlocked<DT_, IT_, 2> zero1;
    DenseVectorBlocked<DT_, IT_, 2> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    DenseVectorBlocked<DT_, IT_, 2> a(10, DT_(7));
    TEST_CHECK_EQUAL(a, a);
    TEST_CHECK_EQUAL(a.bytes(), 20 * sizeof(DT_) + 1 * sizeof(Index));
    DenseVectorBlocked<DT_, IT_, 2> b(10, DT_(5));
    Tiny::Vector<DT_, 2> tv(42);
    b(7, tv);

    DenseVectorBlocked<DT_, IT_, 2> b_r(b, 4, 6);
    TEST_CHECK_EQUAL(b_r(0)[0], b(0+6)[0]);
    TEST_CHECK_EQUAL(b_r(0)[1], b(0+6)[1]);
    TEST_CHECK_EQUAL(b_r(1)[0], b(1+6)[0]);
    TEST_CHECK_EQUAL(b_r(1)[1], b(1+6)[1]);

    DenseVectorBlocked<DT_, IT_, 2> c(b.clone());
    TEST_CHECK_EQUAL(c.size(), b.size());
    //TEST_CHECK_EQUAL(c(7), b(7));
    auto t1 = c(7);
    auto t2 = b(7);
    auto tp1 = b.elements();
    for (Index i(0) ; i < 2 ; ++i)
    {
      TEST_CHECK_EQUAL(t1.v[i], t2.v[i]);
      TEST_CHECK_EQUAL(tp1[7].v[i], t2.v[i]);
    }
    TEST_CHECK_EQUAL(c, b);
    c.convert(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    //TEST_CHECK_EQUAL(c(7), b(7));
    t1 = c(7);
    t2 = b(7);
    for (Index i(0) ; i < 2 ; ++i)
      TEST_CHECK_EQUAL(t1.v[i], t2.v[i]);
    TEST_CHECK_EQUAL(c, b);
    DenseVectorBlocked<float, unsigned int, 2> d;
    d.convert(c);
    DenseVectorBlocked<float, unsigned int, 2> e;
    e.convert(b);
    TEST_CHECK_EQUAL(e.size(), d.size());
    //TEST_CHECK_EQUAL(e(7), d(7));
    t1 = c(7);
    t2 = b(7);
    for (Index i(0) ; i < 2 ; ++i)
      TEST_CHECK_EQUAL(t1.v[i], t2.v[i]);
    TEST_CHECK_EQUAL(e, d);

    b.clone(a);
    TEST_CHECK_NOT_EQUAL((void*)b.elements(), (void*)a.elements());
    c.convert(a);
    TEST_CHECK_EQUAL((void*)c.elements(), (void*)a.elements());
    TEST_CHECK_EQUAL(b, c);
    Tiny::Vector<DT_, 2> tv2(23);
    a(3, tv2);
    TEST_CHECK_EQUAL(a, c);
    TEST_CHECK_NOT_EQUAL(a, b);

    DenseVector<DT_, IT_> dv(12, DT_(2));
    dv(7, DT_(3));
    DenseVectorBlocked<DT_, IT_, 3> f(dv);
    auto t3 = f(2);
    TEST_CHECK_EQUAL(t3.v[0], DT_(2));
    TEST_CHECK_EQUAL(t3.v[1], DT_(3));
    TEST_CHECK_EQUAL((void*)f.elements(), (void*)dv.elements());
    DenseVector<DT_, IT_> dv2(f);
    TEST_CHECK_EQUAL(dv2, dv);
    TEST_CHECK_EQUAL((void*)dv2.elements(), (void*)dv.elements());

    DenseVectorBlocked<DT_, IT_, 3> g(f.size(), f.template elements<Perspective::pod>());
    TEST_CHECK_EQUAL(g, f);
    TEST_CHECK_EQUAL((void*)g.template elements<Perspective::pod>(), (void*)f.template elements<Perspective::pod>());

    // random constructor check
    Random::SeedType seed(Random::SeedType(time(nullptr)));
    std::cout << "seed: " << seed << std::endl;
    Random rng(seed);
    DT_ rnd_range[2];
    IT_ rnd_size = 3*1234;
    rnd_range[0] = DT_(-10);
    rnd_range[1] = DT_(+10);
    DenseVectorBlocked<DT_, IT_, 3> rnd_vec(rng, rnd_size, rnd_range[0], rnd_range[1]);
    TEST_CHECK_EQUAL(rnd_vec.size(), rnd_size);
    DT_ rnd_max = rnd_vec.max_abs_element();
    TEST_CHECK_IN_RANGE(rnd_max, rnd_range[0], rnd_range[1]);
    rnd_vec.scale(rnd_vec, DT_(-1));
    DT_ rnd_min = -rnd_vec.max_abs_element();
    TEST_CHECK_IN_RANGE(rnd_min, rnd_range[0], rnd_range[1]);

  }
};
DenseVectorBlockedTest <float, std::uint32_t> cpu_dv_blocked_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedTest <double, std::uint32_t> cpu_dv_blocked_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedTest <float, std::uint64_t> cpu_dv_blocked_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedTest <double, std::uint64_t> cpu_dv_blocked_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedTest <float, std::uint64_t> mkl_cpu_dv_blocked_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedTest <double, std::uint64_t> mkl_cpu_dv_blocked_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedTest <__float128, std::uint32_t> cpu_dv_blocked_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedTest <__float128, std::uint64_t> cpu_dv_blocked_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedTest <Half, std::uint32_t> dv_blocked_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedTest <Half, std::uint64_t> dv_blocked_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedTest <float, std::uint32_t> cuda_dv_blocked_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedTest <double, std::uint32_t> cuda_dv_blocked_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedTest <float, std::uint64_t> cuda_dv_blocked_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedTest <double, std::uint64_t> cuda_dv_blocked_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorBlockedSerializeTest
  : public UnitTest
{
public:
  DenseVectorBlockedSerializeTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedSerializeTest()
  {
  }

  virtual void run() const override
  {
    DenseVector<DT_, IT_> dv(12, DT_(2));
    dv(7, DT_(3));
    DenseVectorBlocked<DT_, IT_, 3> g(dv);

    std::stringstream mts;
    g.write_out(FileMode::fm_mtx, mts);
    DenseVectorBlocked<DT_, IT_, 3> l(FileMode::fm_mtx, mts);
    TEST_CHECK_EQUAL(l, g);
    //for (Index i(0) ; i < g.template size<Perspective::pod>() ; ++i)
    //  TEST_CHECK_EQUAL_WITHIN_EPS(l.template elements<Perspective::pod>()[i], g.template elements<Perspective::pod>()[i], DT_(1e-4));

    std::stringstream ts;
    g.write_out(FileMode::fm_exp, ts);
    DenseVectorBlocked<DT_, IT_, 3> m(FileMode::fm_exp, ts);
    TEST_CHECK_EQUAL(m, g);
    //for (Index i(0) ; i < k.size() ; ++i)
    //  TEST_CHECK_EQUAL_WITHIN_EPS(m(i), k(i), DT_(1e-4));

    BinaryStream bs;
    g.write_out(FileMode::fm_dvb, bs);
    bs.seekg(0);
    DenseVectorBlocked<DT_, IT_, 3> n(FileMode::fm_dvb, bs);
    TEST_CHECK_EQUAL(n, g);
    //for (Index i(0) ; i < k.size() ; ++i)
    //  TEST_CHECK_EQUAL_WITHIN_EPS(n(i), k(i), DT_(1e-5));

    DenseVector<DT_,IT_> zfp_tester(20);
    for(Index i(0); i < zfp_tester.size() ; ++i)
    {
      zfp_tester(i, DT_(7)/(((DT_(i)+DT_(0.5))*DT_(i+1))));
    }
    DenseVectorBlocked<DT_, IT_, 5> tester_blocked(zfp_tester);
    auto op = g.serialize(LAFEM::SerialConfig(false, false));
    DenseVectorBlocked<DT_, IT_, 3> o(op);
    TEST_CHECK_EQUAL(o, g);
#ifdef FEAT_HAVE_ZLIB
    auto zb = tester_blocked.serialize(LAFEM::SerialConfig(true,false));
    DenseVectorBlocked<DT_, IT_, 5> zlib(zb);
    TEST_CHECK_EQUAL(zlib, tester_blocked);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zp = tester_blocked.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    DenseVectorBlocked<DT_, IT_, 5> zfp(zp);
    for (Index i(0) ; i < tester_blocked.size() ; ++i)
    {
      auto tst = zfp(i);
      for(int j(0) ; j < tst.n ; ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(tst[j], tester_blocked(i)[j], DT_(1e-5));
      }
    }
#endif
  }
};
DenseVectorBlockedSerializeTest <float, std::uint32_t> cpu_dense_vector_blocked_serialize_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedSerializeTest <double, std::uint32_t> cpu_dense_vector_blocked_serialize_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedSerializeTest <float, std::uint64_t> cpu_dense_vector_blocked_serialize_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedSerializeTest <double, std::uint64_t> cpu_dense_vector_blocked_serialize_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedSerializeTest <float, std::uint64_t> mkl_cpu_dense_vector_blocked_serialize_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedSerializeTest <double, std::uint64_t> mkl_cpu_dense_vector_blocked_serialize_test_double_uint64(PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//DenseVectorBlockedSerializeTest <__float128, std::uint32_t> cpu_dense_vector_blocked_serialize_test_float128_uint32(PreferredBackend::generic);
//DenseVectorBlockedSerializeTest <__float128, std::uint64_t> cpu_dense_vector_blocked_serialize_test_float128_uint64(PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedSerializeTest <Half, std::uint32_t> dv_blocked_serialize_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedSerializeTest <Half, std::uint64_t> dv_blocked_serialize_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedSerializeTest <float, std::uint32_t> cuda_dense_vector_blocked_serialize_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedSerializeTest <double, std::uint32_t> cuda_dense_vector_blocked_serialize_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedSerializeTest <float, std::uint64_t> cuda_dense_vector_blocked_serialize_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedSerializeTest <double, std::uint64_t> cuda_dense_vector_blocked_serialize_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedAxpyTest
  : public UnitTest
{
public:
  DenseVectorBlockedAxpyTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedAxpyTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    DT_ s(DT_(0.4711));
    Tiny::Vector<DT_, block_size_> bs;
    for(int j(0); j < block_size_; ++j)
      bs[j] = s + DT_(j);

    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref_bs(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, block_size_> tv1;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv1.v[j] = DT_(i % 100) * DT_(1.234) + DT_(j);
        a(i, tv1);
        Tiny::Vector<DT_, block_size_> tv2;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv2.v[j] = DT_(2) - DT_(i % (42 + j));
        b(i, tv2);
        ref(i, s * a(i) + b(i));

        ref_bs(i, Tiny::component_product(bs, a(i)) + ref(i));
      }

      DenseVectorBlocked<DT_, IT_, block_size_> c(size);
      c.axpy(a, b, s);
      for (Index i(0) ; i < size ; ++i)
        for (Index j(0) ; j < block_size_ ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(c(i).v[j], ref(i).v[j], DT_(eps));

      DenseVectorBlocked<DT_, IT_, block_size_> a_tmp(size);
      a_tmp.copy(a);
      a.axpy(a, b, s);
      for (Index i(0) ; i < size ; ++i)
        for (Index j(0) ; j < block_size_ ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(a(i).v[j], ref(i).v[j], DT_(eps));

      a.copy(a_tmp);
      b.axpy(a, b, s);
      for (Index i(0) ; i < size ; ++i)
        for (Index j(0) ; j < block_size_ ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(b(i).v[j], ref(i).v[j], DT_(eps));

      a.copy(a_tmp);
      b.axpy_blocked(a, b, bs);
      for (Index i(0) ; i < size ; ++i)
        for (Index j(0) ; j < block_size_ ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(b(i).v[j], ref_bs(i).v[j], DT_(eps));
    }
  }
};
DenseVectorBlockedAxpyTest <float, std::uint32_t, 2> dv_axpy_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedAxpyTest <double, std::uint32_t, 2> dv_axpy_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedAxpyTest <float, std::uint64_t, 3> dv_axpy_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedAxpyTest <double, std::uint64_t, 3> dv_axpy_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedAxpyTest <float, std::uint64_t, 2> mkl_cpu_dense_vector_blocked_axpy_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedAxpyTest <double, std::uint64_t, 3> mkl_cpu_dense_vector_blocked_axpy_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedAxpyTest <__float128, std::uint32_t, 2> dv_axpy_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedAxpyTest <__float128, std::uint64_t, 3> dv_axpy_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedAxpyTest <Half, std::uint32_t, 2> dv_blocked_axpy_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedAxpyTest <Half, std::uint64_t, 3> dv_blocked_axpy_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedAxpyTest <float, std::uint32_t, 2> cuda_dv_axpy_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedAxpyTest <double, std::uint32_t, 2> cuda_dv_axpy_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedAxpyTest <float, std::uint64_t, 3> cuda_dv_axpy_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedAxpyTest <double, std::uint64_t, 3> cuda_dv_axpy_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedDotTest
  : public UnitTest
{
public:
  DenseVectorBlockedDotTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedDotTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedDotTest()
  {
  }

  virtual void run() const override
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.7));

    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> a2(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      const DT_ den(DT_(1) / DT_(size * block_size_));
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, block_size_> tv1;
        Tiny::Vector<DT_, block_size_> tv3;
        for (Index j(0) ; j < block_size_ ; ++j)
        {
          tv1.v[j] = DT_(i * block_size_ + j + 1) * den;
          tv3.v[j] = Math::sqrt(DT_(i))*j;
        }
        a(i, tv1);
        a2(i, tv3);

        Tiny::Vector<DT_, block_size_> tv2;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv2.v[j] = DT_(1) / DT_(i * block_size_ + j + 1);
        b(i, tv2);
      }

      Tiny::Vector<DT_, block_size_> d(a.dot_blocked(b));
      Tiny::Vector<DT_, block_size_> e(b.dot_blocked(a));
      Tiny::Vector<DT_, block_size_> f(a2.dot_blocked(a2));

      Tiny::Vector<DT_, block_size_> ref_bs;
      Tiny::Vector<DT_, block_size_> ref2_bs;
      for (Index j(0) ; j < block_size_ ; ++j)
      {
        ref_bs.v[j] = DT_(size)/DT_(size*block_size_);
        ref2_bs.v[j] = DT_(j*j)*(DT_(size*size)-DT_(size))/DT_(2);
      }

      for (Index j(0) ; j < block_size_ ; ++j)
      {
       TEST_CHECK_EQUAL_WITHIN_EPS(ref_bs.v[j], d.v[j], eps);
       TEST_CHECK_EQUAL_WITHIN_EPS(ref_bs.v[j], e.v[j], eps);
       TEST_CHECK_EQUAL_WITHIN_EPS(ref2_bs.v[j], f.v[j], eps);
      }

      // a*b = 1
      DT_ ref(DT_(1));
      DT_ c  = a.dot(b);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = b.dot(a);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = b.dot(b);
      ref = b.norm2();
      ref *= ref;
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);

    }
  }
};
DenseVectorBlockedDotTest <float, std::uint32_t, 3> dv_dot_product_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedDotTest <double, std::uint32_t, 3> dv_dot_product_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedDotTest <float, std::uint64_t, 2> dv_dot_product_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedDotTest <double, std::uint64_t, 2> dv_dot_product_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedDotTest <float, std::uint64_t, 3> mkl_cpu_dense_vector_blocked_dot_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedDotTest <double, std::uint64_t, 2> mkl_cpu_dense_vector_blocked_dot_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedDotTest <__float128, std::uint32_t, 3> dv_dot_product_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedDotTest <__float128, std::uint64_t, 2> dv_dot_product_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedDotTest <Half, std::uint32_t,3 > dv_blocked_dot_product_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedDotTest <Half, std::uint64_t, 2> dv_blocked_dot_product_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedDotTest <float, std::uint32_t, 3> cuda_dv_dot_product_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedDotTest <double, std::uint32_t, 3> cuda_dv_dot_product_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedDotTest <float, std::uint64_t, 2> cuda_dv_dot_product_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedDotTest <double, std::uint64_t, 2> cuda_dv_dot_product_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedTripleDotTest
  : public UnitTest
{
public:
  DenseVectorBlockedTripleDotTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedTripleDotTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedTripleDotTest()
  {
  }

  virtual void run() const override
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      DenseVectorBlocked<DT_, IT_, block_size_> c(size);
      Tiny::Vector<DT_, block_size_> ref_bs(0);
      const DT_ den(DT_(1) / DT_(size * block_size_));
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, block_size_> tv1;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv1.v[j] = DT_(i * block_size_ + j + 1) * den;
        a(i, tv1);
        Tiny::Vector<DT_, block_size_> tv2;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv2.v[j] = DT_(1) / DT_(i * block_size_ + j + 1);
        b(i, tv2);
        Tiny::Vector<DT_, block_size_> tv3;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv3.v[j] =  DT_(i * block_size_ + j + 1);
        c(i, tv3);
      }

      for(int j(0); j<block_size_; ++j)
        ref_bs.v[j] = DT_(size-1)/DT_(2)+DT_(j+1)/DT_(block_size_);

      Tiny::Vector<DT_, block_size_> d(a.triple_dot_blocked(b, c));
      Tiny::Vector<DT_, block_size_> e(b.triple_dot_blocked(a, c));
      Tiny::Vector<DT_, block_size_> f(c.triple_dot_blocked(a, b));

      for(int j(0); j<block_size_; ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(d.v[j], ref_bs.v[j], eps*DT_(250));
        TEST_CHECK_EQUAL_WITHIN_EPS(e.v[j], ref_bs.v[j], eps*DT_(250));
        TEST_CHECK_EQUAL_WITHIN_EPS(f.v[j], ref_bs.v[j], eps*DT_(250));
      }

      DenseVector<DT_, IT_> ref_a;
      ref_a.convert(a);
      DenseVector<DT_, IT_> ref_b;
      ref_b.convert(b);
      DenseVector<DT_, IT_> ref_c;
      ref_c.convert(c);

      DT_ ref(ref_a.triple_dot(ref_b, ref_c));
      DT_ res  = a.triple_dot(b, c);
      TEST_CHECK_EQUAL_WITHIN_EPS(res, ref, eps);
      res = b.triple_dot(a, c);
      TEST_CHECK_EQUAL_WITHIN_EPS(res, ref, eps);
      res = c.triple_dot(b, a);
      TEST_CHECK_EQUAL_WITHIN_EPS(res, ref, eps);
    }
  }
};
DenseVectorBlockedTripleDotTest <float, std::uint32_t, 3> dv_triple_dot_product_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedTripleDotTest <double, std::uint32_t, 3> dv_triple_dot_product_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedTripleDotTest <float, std::uint64_t, 2> dv_triple_dot_product_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedTripleDotTest <double, std::uint64_t, 2> dv_triple_dot_product_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedTripleDotTest <float, std::uint64_t,3 > mkl_dv_blocked_triple_dot_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedTripleDotTest <double, std::uint64_t, 2> mkl_dv_blocked_triple_dot_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedTripleDotTest <__float128, std::uint32_t, 3> dv_triple_dot_product_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedTripleDotTest <__float128, std::uint64_t, 2> dv_triple_dot_product_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedTripleDotTest <Half, std::uint32_t, 3> dv_blocked_triple_dot_product_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedTripleDotTest <Half, std::uint64_t, 2> dv_blocked_triple_dot_product_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedTripleDotTest <float, std::uint32_t, 3> cuda_dv_triple_dot_product_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedTripleDotTest <double, std::uint32_t, 3> cuda_dv_triple_dot_product_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedTripleDotTest <float, std::uint64_t, 2> cuda_dv_triple_dot_product_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedTripleDotTest <double, std::uint64_t, 2> cuda_dv_triple_dot_product_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedComponentProductTest
  : public UnitTest
{
public:
  DenseVectorBlockedComponentProductTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedComponentProductTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedComponentProductTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.7));

      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref2(size);

      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, block_size_> tv1;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv1.v[j] = DT_((DT_(j) + DT_(i)) * DT_(1.234) / DT_(size));
        a(i, tv1);
        Tiny::Vector<DT_, block_size_> tv2;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv2.v[j] = DT_(size*2*block_size_ - i + 2*j);
        b(i, tv2);
        ref(i, component_product(a(i), b(i)) );
        ref2(i, component_product(a(i), a(i)) );
      }

      DenseVectorBlocked<DT_, IT_, block_size_> a_tmp(size);
      a_tmp.copy(a);
      DenseVectorBlocked<DT_, IT_, block_size_> b_tmp(size);
      b_tmp.copy(b);
      DenseVectorBlocked<DT_, IT_, block_size_> c(size);

      c.component_product(a, b);
      for (Index i(0); i < c.template size<Perspective::pod>(); ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(c.template elements<Perspective::pod>()[i], ref.template elements<Perspective::pod>()[i], eps);
      }

      a.component_product(a, b);
      for (Index i(0); i < a.template size<Perspective::pod>(); ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(a.template elements<Perspective::pod>()[i], ref.template elements<Perspective::pod>()[i], eps);
      }

      a.copy(a_tmp);
      b.component_product(a, b);
      for (Index i(0); i < a.template size<Perspective::pod>(); ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(b.template elements<Perspective::pod>()[i], ref.template elements<Perspective::pod>()[i], eps);
      }

      a.component_product(a, a);
      for (Index i(0); i < a.template size<Perspective::pod>(); ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(a.template elements<Perspective::pod>()[i], ref2.template elements<Perspective::pod>()[i], eps);
      }
    }
  }
};
DenseVectorBlockedComponentProductTest <float, std::uint32_t, 3> dv_blocked_component_product_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedComponentProductTest <double, std::uint32_t, 3> dv_blocked_component_product_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedComponentProductTest <float, std::uint64_t, 2> dv_blocked_component_product_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedComponentProductTest <double, std::uint64_t, 2> dv_blocked_component_product_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedComponentProductTest <float, std::uint64_t, 3> mkl_dv_blocked_component_product_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedComponentProductTest <double, std::uint64_t, 2> mkl_dv_blocked_component_product_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedComponentProductTest <__float128, std::uint32_t, 3> dv_blocked_component_product_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedComponentProductTest <__float128, std::uint64_t, 2> dv_blocked_component_product_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedComponentProductTest <Half, std::uint32_t, 3> dv_blocked_component_product_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedComponentProductTest <Half, std::uint64_t, 2> dv_blocked_component_product_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedComponentProductTest <float, std::uint32_t, 3> cuda_dv_blocked_component_product_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedComponentProductTest <double, std::uint32_t, 3> cuda_dv_blocked_component_product_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedComponentProductTest <float, std::uint64_t, 2> cuda_dv_blocked_component_product_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedComponentProductTest <double, std::uint64_t, 2> cuda_dv_blocked_component_product_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedScaleTest
  : public UnitTest
{
public:
  DenseVectorBlockedScaleTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedScaleTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedScaleTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DT_ s(DT_(4.321));
      Tiny::Vector<DT_, block_size_> bs;
      for(Index j(0); j < block_size_; ++j)
        bs[j] = s + DT_(j);
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref_bs(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, block_size_> tv1;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv1.v[j] = DT_(DT_(i) * DT_(block_size_) + DT_(j) * DT_(1.234));
        a(i, tv1);
        ref(i, a(i) * s);

        ref_bs(i, Tiny::component_product(bs, s*a(i)));
      }

      DenseVectorBlocked<DT_, IT_, block_size_> b(size);

      b.scale(a, s);
      TEST_CHECK_EQUAL(b, ref);

      a.scale(a, s);
      TEST_CHECK_EQUAL(a, ref);

      a.scale_blocked(a, bs);
      TEST_CHECK_EQUAL(a, ref_bs);
    }
  }
};
DenseVectorBlockedScaleTest <float, std::uint32_t, 3> dv_scale_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedScaleTest <double, std::uint32_t, 3> dv_scale_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedScaleTest <float, std::uint64_t, 2> dv_scale_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedScaleTest <double, std::uint64_t, 2> dv_scale_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedScaleTest <float, std::uint64_t, 3> mkl_dv_blocked_scale_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedScaleTest <double, std::uint64_t, 2> mkl_dv_blocked_scale_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedScaleTest <__float128, std::uint32_t, 3> dv_scale_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedScaleTest <__float128, std::uint64_t, 2> dv_scale_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedScaleTest <Half, std::uint32_t, 3> dv_blocked_scale_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedScaleTest <Half, std::uint64_t, 2> dv_blocked_scale_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedScaleTest <float, std::uint32_t, 3> cuda_dv_scale_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedScaleTest <double, std::uint32_t, 3> cuda_dv_scale_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedScaleTest <float, std::uint64_t, 2> cuda_dv_scale_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedScaleTest <double, std::uint64_t, 2> cuda_dv_scale_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedNorm2Test
  : public UnitTest
{
public:
  DenseVectorBlockedNorm2Test(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedNorm2Test", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedNorm2Test()
  {
  }

  virtual void run() const override
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      Tiny::Vector<DT_, block_size_> ref_bs(0);
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);

      for (Index i(0) ; i < size ; ++i)
      {
        // a[i] = 1/sqrt(2^i) = (1/2)^(i/2)
        Tiny::Vector<DT_, block_size_> tv1;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv1.v[j] = Math::pow(DT_(0.5), DT_(0.5) * DT_(i * block_size_ + j));
        a(i, tv1);
      }

      // ||a||_2 = sqrt(2 - 2^{1-n})
      const DT_ ref(Math::sqrt(DT_(2) - Math::pow(DT_(0.5), DT_(size*block_size_-1))));
      for (Index j(0) ; j < block_size_ ; ++j)
      {
        ref_bs.v[j] = Math::pow(DT_(0.5), DT_(j))*(DT_(1)-Math::pow(DT_(0.5), DT_(block_size_*size)))/(DT_(1)-Math::pow(DT_(0.5), DT_(block_size_)));
        ref_bs.v[j] = Math::sqrt(ref_bs.v[j]);
      }

      DT_ c = a.norm2();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);

      c = a.norm2sqr();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref*ref, eps);

      Tiny::Vector<DT_, block_size_> d(a.norm2_blocked());
      for(Index j(0); j < block_size_; ++j)
        TEST_CHECK_EQUAL_WITHIN_EPS(d.v[j], ref_bs.v[j], eps);

      Tiny::Vector<DT_, block_size_> e(a.norm2sqr_blocked());
      for(Index j(0); j < block_size_; ++j)
        TEST_CHECK_EQUAL_WITHIN_EPS(e.v[j], ref_bs.v[j]*ref_bs.v[j], eps);
    }
  }
};
DenseVectorBlockedNorm2Test <float, std::uint32_t, 2> dv_norm2_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedNorm2Test <double, std::uint32_t, 2> dv_norm2_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedNorm2Test <float, std::uint64_t, 3> dv_norm2_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedNorm2Test <double, std::uint64_t, 3> dv_norm2_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedNorm2Test <float, std::uint64_t, 2> mkl_dv_norm2_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedNorm2Test <double, std::uint64_t, 3> mkl_dv_norm2_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedNorm2Test <__float128, std::uint32_t, 2> dv_norm2_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedNorm2Test <__float128, std::uint64_t, 3> dv_norm2_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedNorm2Test <Half, std::uint32_t, 2> dv_blocked_norm2_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedNorm2Test <Half, std::uint64_t, 3> dv_blocked_norm2_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedNorm2Test <float, std::uint32_t, 2> cuda_dv_norm2_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedNorm2Test <double, std::uint32_t, 2> cuda_dv_norm2_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedNorm2Test <float, std::uint64_t, 3> cuda_dv_norm2_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedNorm2Test <double, std::uint64_t, 3> cuda_dv_norm2_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMaxAbsElementTest
  : public UnitTest
{
public:
  DenseVectorBlockedMaxAbsElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMaxAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedMaxAbsElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e4) ; size*=2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, block_size_> tv1;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv1.v[j]  = DT_(DT_((i * block_size_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        a(i, tv1);
      }

      Tiny::Vector<DT_, block_size_> b(a.max_abs_element_blocked());
      Tiny::Vector<DT_, block_size_> ref_bs;
      for (Index j(0) ; j < block_size_ ; ++j)
        ref_bs.v[j] = DT_((size-1) * block_size_  + j);

      for(Index j(0); j < block_size_; ++j)
        TEST_CHECK_EQUAL(b.v[j], ref_bs.v[j]);

      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);
      DT_ max = a.max_abs_element();

      TEST_CHECK_EQUAL(max, DT_((size*block_size_) -1));
    }
  }
};
DenseVectorBlockedMaxAbsElementTest <float, std::uint32_t, 2> dv_max_abs_element_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedMaxAbsElementTest <double, std::uint32_t, 2> dv_max_abs_element_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedMaxAbsElementTest <float, std::uint64_t, 3> dv_max_abs_element_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedMaxAbsElementTest <double, std::uint64_t, 3> dv_max_abs_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedMaxAbsElementTest <float, std::uint64_t, 2> mkl_dv_max_abs_element_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedMaxAbsElementTest <double, std::uint64_t, 3> mkl_dv_max_abs_element_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedMaxAbsElementTest <__float128, std::uint32_t, 2> dv_max_abs_element_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedMaxAbsElementTest <__float128, std::uint64_t, 3> dv_max_abs_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedMaxAbsElementTest <Half, std::uint32_t, 2> dv_blocked_max_abs_element_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedMaxAbsElementTest <Half, std::uint64_t, 3> dv_blocked_max_abs_element_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedMaxAbsElementTest <float, std::uint32_t, 2> cuda_dv_max_abs_element_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedMaxAbsElementTest <double, std::uint32_t, 2> cuda_dv_max_abs_element_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedMaxAbsElementTest <float, std::uint64_t, 3> cuda_dv_max_abs_element_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedMaxAbsElementTest <double, std::uint64_t, 3> cuda_dv_max_abs_element_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMinAbsElementTest
  : public UnitTest
{
public:
  DenseVectorBlockedMinAbsElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMinAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedMinAbsElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e4) ; size*=2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, block_size_> tv1;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv1.v[j]  = DT_(DT_((i * block_size_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        a(i, tv1);
      }

      Tiny::Vector<DT_, block_size_> b(a.min_abs_element_blocked());
      Tiny::Vector<DT_, block_size_> ref_bs;
      for (Index j(0) ; j < block_size_ ; ++j)
        ref_bs.v[j] = DT_(j);

      for(Index j(0); j < block_size_; ++j)
        TEST_CHECK_EQUAL(b.v[j], ref_bs.v[j]);

      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ min = a.min_abs_element();

      TEST_CHECK_EQUAL(min, DT_(0));
    }
  }
};
DenseVectorBlockedMinAbsElementTest <float, std::uint32_t, 2> dv_min_abs_element_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedMinAbsElementTest <double, std::uint32_t, 2> dv_min_abs_element_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedMinAbsElementTest <float, std::uint64_t, 3> dv_min_abs_element_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedMinAbsElementTest <double, std::uint64_t, 3> dv_min_abs_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedMinAbsElementTest <float, std::uint64_t, 2> mkl_dv_min_abs_element_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedMinAbsElementTest <double, std::uint64_t, 3> mkl_dv_min_abs_element_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedMinAbsElementTest <__float128, std::uint32_t, 2> dv_min_abs_element_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedMinAbsElementTest <__float128, std::uint64_t, 3> dv_min_abs_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedMinAbsElementTest <Half, std::uint32_t, 2> dv_blocked_min_abs_element_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedMinAbsElementTest <Half, std::uint64_t, 3> dv_blocked_min_abs_element_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedMinAbsElementTest <float, std::uint32_t, 2> cuda_dv_min_abs_element_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedMinAbsElementTest <double, std::uint32_t, 2> cuda_dv_min_abs_element_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedMinAbsElementTest <float, std::uint64_t, 3> cuda_dv_min_abs_element_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedMinAbsElementTest <double, std::uint64_t, 3> cuda_dv_min_abs_element_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMaxElementTest
  : public UnitTest
{
public:
  DenseVectorBlockedMaxElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMaxElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedMaxElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e4) ; size*=2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, block_size_> tv1;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv1.v[j]  = DT_(DT_((i * block_size_ + j)));
        a(i, tv1);
      }

      Tiny::Vector<DT_, block_size_> b(a.max_element_blocked());
      Tiny::Vector<DT_, block_size_> ref_bs;
      for (Index j(0) ; j < block_size_ ; ++j)
        ref_bs.v[j] = DT_((size-1)*block_size_+j);

      for(Index j(0); j < block_size_; ++j)
        TEST_CHECK_EQUAL(b.v[j], ref_bs.v[j]);

      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ max = a.max_element();

      TEST_CHECK_EQUAL(max, DT_((size*block_size_) -1));
    }
  }
};
DenseVectorBlockedMaxElementTest <float, std::uint32_t, 2> dv_max_element_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedMaxElementTest <double, std::uint32_t, 2> dv_max_element_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedMaxElementTest <float, std::uint64_t, 3> dv_max_element_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedMaxElementTest <double, std::uint64_t, 3> dv_max_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedMaxElementTest <float, std::uint64_t, 2> mkl_dv_max_element_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedMaxElementTest <double, std::uint64_t, 3> mkl_dv_max_element_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedMaxElementTest <__float128, std::uint32_t, 2> dv_max_element_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedMaxElementTest <__float128, std::uint64_t, 3> dv_max_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedMaxElementTest <Half, std::uint32_t, 2> dv_blocked_max_element_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedMaxElementTest <Half, std::uint64_t, 3> dv_blocked_max_element_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedMaxElementTest <float, std::uint32_t, 2> cuda_dv_max_element_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedMaxElementTest <double, std::uint32_t, 2> cuda_dv_max_element_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedMaxElementTest <float, std::uint64_t, 3> cuda_dv_max_element_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedMaxElementTest <double, std::uint64_t, 3> cuda_dv_max_element_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMinElementTest
  : public UnitTest
{
public:
  DenseVectorBlockedMinElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMinElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedMinElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e4) ; size*=2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, block_size_> tv1;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv1.v[j]  = DT_(DT_((i * block_size_ + j)));
        a(i, tv1);
      }

      Tiny::Vector<DT_, block_size_> b(a.min_element_blocked());
      Tiny::Vector<DT_, block_size_> ref_bs;
      for (Index j(0) ; j < block_size_ ; ++j)
        ref_bs.v[j] = DT_(j);

      for(Index j(0); j < block_size_; ++j)
        TEST_CHECK_EQUAL(b.v[j], ref_bs.v[j]);

      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ min = a.min_element();

      TEST_CHECK_EQUAL(min, DT_(0));
    }
  }
};
DenseVectorBlockedMinElementTest <float, std::uint32_t, 2> dv_min_element_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedMinElementTest <double, std::uint32_t, 2> dv_min_element_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedMinElementTest <float, std::uint64_t, 3> dv_min_element_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedMinElementTest <double, std::uint64_t, 3> dv_min_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedMinElementTest <float, std::uint64_t, 2> mkl_dv_min_element_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedMinElementTest <double, std::uint64_t, 3> mkl_dv_min_element_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedMinElementTest <__float128, std::uint32_t, 2> dv_min_element_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedMinElementTest <__float128, std::uint64_t, 3> dv_min_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedMinElementTest <Half, std::uint32_t, 2> dv_blocked_min_element_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedMinElementTest <Half, std::uint64_t, 3> dv_blocked_min_element_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedMinElementTest <float, std::uint32_t, 2> cuda_dv_min_element_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedMinElementTest <double, std::uint32_t, 2> cuda_dv_min_element_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedMinElementTest <float, std::uint64_t, 3> cuda_dv_min_element_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedMinElementTest <double, std::uint64_t, 3> cuda_dv_min_element_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedPermuteTest
  : public UnitTest
{
public:
  DenseVectorBlockedPermuteTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedPermuteTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorBlockedPermuteTest()
  {
  }

  virtual void run() const override
  {
    bool one_error(false);

    for (Index size(10) ; size < Index(1e4) ; size*=2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, block_size_> tv1;
        for (Index j(0) ; j < block_size_ ; ++j)
          tv1.v[j]  = DT_(DT_((i * block_size_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        a(i, tv1);
      }

      DenseVectorBlocked<DT_, IT_, block_size_> ref(a.size());
      ref.copy(a);

      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      //allow identity permutation once per test due to rng...
      if (a==ref)
      {
        if(one_error)
        {
          TEST_CHECK(false);
        }
        else
        {
          one_error = true;
        }
      }

      auto prm_inv = prm_rnd.inverse();
      a.permute(prm_inv);
      TEST_CHECK_EQUAL(a, ref);
    }
  }
};
DenseVectorBlockedPermuteTest <float, std::uint32_t, 2> dv_permute_test_float_uint32(PreferredBackend::generic);
DenseVectorBlockedPermuteTest <double, std::uint32_t, 2> dv_permute_test_double_uint32(PreferredBackend::generic);
DenseVectorBlockedPermuteTest <float, std::uint64_t, 3> dv_permute_test_float_uint64(PreferredBackend::generic);
DenseVectorBlockedPermuteTest <double, std::uint64_t, 3> dv_permute_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseVectorBlockedPermuteTest <float, std::uint64_t, 2> mkl_dv_permute_test_float_uint64(PreferredBackend::mkl);
DenseVectorBlockedPermuteTest <double, std::uint64_t, 3> mkl_dv_permute_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedPermuteTest <__float128, std::uint32_t, 2> dv_permute_test_float128_uint32(PreferredBackend::generic);
DenseVectorBlockedPermuteTest <__float128, std::uint64_t, 3> dv_permute_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorBlockedPermuteTest <Half, std::uint32_t, 2> dv_blocked_permute_test_half_uint32(PreferredBackend::generic);
DenseVectorBlockedPermuteTest <Half, std::uint64_t, 3> dv_blocked_permute_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedPermuteTest <float, std::uint32_t, 2> cuda_dv_permute_test_float_uint32(PreferredBackend::cuda);
DenseVectorBlockedPermuteTest <double, std::uint32_t, 2> cuda_dv_permute_test_double_uint32(PreferredBackend::cuda);
DenseVectorBlockedPermuteTest <float, std::uint64_t, 3> cuda_dv_permute_test_float_uint64(PreferredBackend::cuda);
DenseVectorBlockedPermuteTest <double, std::uint64_t, 3> cuda_dv_permute_test_double_uint64(PreferredBackend::cuda);
#endif
