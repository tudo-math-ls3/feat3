// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
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
class DenseVectorBlockedTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedTest()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedTest")
  {
  }

  virtual ~DenseVectorBlockedTest()
  {
  }

  virtual void run() const override
  {
    DenseVectorBlocked<Mem_, DT_, IT_, 2> zero1;
    DenseVectorBlocked<Mem::Main, DT_, IT_, 2> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    DenseVectorBlocked<Mem_, DT_, IT_, 2> a(10, DT_(7));
    TEST_CHECK_EQUAL(a, a);
    TEST_CHECK_EQUAL(a.bytes(), 20 * sizeof(DT_) + 1 * sizeof(Index));
    DenseVectorBlocked<Mem_, DT_, IT_, 2> b(10, DT_(5));
    Tiny::Vector<DT_, 2> tv(42);
    b(7, tv);

    DenseVectorBlocked<Mem_, DT_, IT_, 2> b_r(b, 4, 6);
    TEST_CHECK_EQUAL(b_r(0)[0], b(0+6)[0]);
    TEST_CHECK_EQUAL(b_r(0)[1], b(0+6)[1]);
    TEST_CHECK_EQUAL(b_r(1)[0], b(1+6)[0]);
    TEST_CHECK_EQUAL(b_r(1)[1], b(1+6)[1]);

    DenseVectorBlocked<Mem_, DT_, IT_, 2> c(b.clone());
    TEST_CHECK_EQUAL(c.size(), b.size());
    //TEST_CHECK_EQUAL(c(7), b(7));
    auto t1 = c(7);
    auto t2 = b(7);
    auto tp1 = b.elements();
    for (Index i(0) ; i < 2 ; ++i)
    {
      TEST_CHECK_EQUAL(t1.v[i], t2.v[i]);
      if (std::is_same<Mem_, Mem::Main>::value)
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
    DenseVectorBlocked<Mem::Main, float, unsigned int, 2> d;
    d.convert(c);
    DenseVectorBlocked<Mem::Main, float, unsigned int, 2> e;
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

    DenseVector<Mem_, DT_, IT_> dv(12, DT_(2));
    dv(7, DT_(3));
    DenseVectorBlocked<Mem_, DT_, IT_, 3> f(dv);
    auto t3 = f(2);
    TEST_CHECK_EQUAL(t3.v[0], DT_(2));
    TEST_CHECK_EQUAL(t3.v[1], DT_(3));
    TEST_CHECK_EQUAL((void*)f.elements(), (void*)dv.elements());
    DenseVector<Mem_, DT_, IT_> dv2(f);
    TEST_CHECK_EQUAL(dv2, dv);
    TEST_CHECK_EQUAL((void*)dv2.elements(), (void*)dv.elements());

    DenseVectorBlocked<Mem_, DT_, IT_, 3> g(f.size(), f.template elements<Perspective::pod>());
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
    DenseVectorBlocked<Mem_, DT_, IT_, 3> rnd_vec(rng, rnd_size, rnd_range[0], rnd_range[1]);
    TEST_CHECK_EQUAL(rnd_vec.size(), rnd_size);
    DT_ rnd_max = rnd_vec.max_abs_element();
    TEST_CHECK_IN_RANGE(rnd_max, rnd_range[0], rnd_range[1]);
    rnd_vec.scale(rnd_vec, DT_(-1));
    DT_ rnd_min = -rnd_vec.max_abs_element();
    TEST_CHECK_IN_RANGE(rnd_min, rnd_range[0], rnd_range[1]);

  }
};
DenseVectorBlockedTest<Mem::Main, float, unsigned int> cpu_dense_vector_blocked_test_float_uint;
DenseVectorBlockedTest<Mem::Main, double, unsigned int> cpu_dense_vector_blocked_test_double_uint;
DenseVectorBlockedTest<Mem::Main, float, unsigned long> cpu_dense_vector_blocked_test_float_ulong;
DenseVectorBlockedTest<Mem::Main, double, unsigned long> cpu_dense_vector_blocked_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedTest<Mem::Main, __float128, unsigned int> cpu_dense_vector_blocked_test_float128_uint;
DenseVectorBlockedTest<Mem::Main, __float128, unsigned long> cpu_dense_vector_blocked_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedTest<Mem::CUDA, float, unsigned int> cuda_dense_vector_blocked_test_float_uint;
DenseVectorBlockedTest<Mem::CUDA, double, unsigned int> cuda_dense_vector_blocked_test_double_uint;
DenseVectorBlockedTest<Mem::CUDA, float, unsigned long> cuda_dense_vector_blocked_test_float_ulong;
DenseVectorBlockedTest<Mem::CUDA, double, unsigned long> cuda_dense_vector_blocked_test_double_ulong;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class DenseVectorBlockedSerializeTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedSerializeTest()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedSerializeTest")
  {
  }

  virtual ~DenseVectorBlockedSerializeTest()
  {
  }

  virtual void run() const override
  {
    DenseVector<Mem_, DT_, IT_> dv(12, DT_(2));
    dv(7, DT_(3));
    DenseVectorBlocked<Mem_, DT_, IT_, 3> g(dv);

    std::stringstream mts;
    g.write_out(FileMode::fm_mtx, mts);
    DenseVectorBlocked<Mem_, DT_, IT_, 3> l(FileMode::fm_mtx, mts);
    TEST_CHECK_EQUAL(l, g);
    /*for (Index i(0) ; i < g.template size<Perspective::pod>() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(l.template elements<Perspective::pod>()[i], g.template elements<Perspective::pod>()[i], DT_(1e-4));*/

    std::stringstream ts;
    g.write_out(FileMode::fm_exp, ts);
    DenseVectorBlocked<Mem_, DT_, IT_, 3> m(FileMode::fm_exp, ts);
    TEST_CHECK_EQUAL(m, g);
    /*for (Index i(0) ; i < k.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(m(i), k(i), DT_(1e-4));*/

    BinaryStream bs;
    g.write_out(FileMode::fm_dvb, bs);
    bs.seekg(0);
    DenseVectorBlocked<Mem_, DT_, IT_, 3> n(FileMode::fm_dvb, bs);
    TEST_CHECK_EQUAL(n, g);
    /*for (Index i(0) ; i < k.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(n(i), k(i), DT_(1e-5));*/

    DenseVector<Mem_,DT_,IT_> zfp_tester(20);
    for(Index i(0); i < zfp_tester.size() ; ++i)
    {
      zfp_tester(i, DT_(7)/(((DT_(i)+DT_(0.5))*DT_(i+1))));
    }
    DenseVectorBlocked<Mem_, DT_, IT_, 5> tester_blocked(zfp_tester);
    auto op = g.serialize(LAFEM::SerialConfig(false, false));
    DenseVectorBlocked<Mem_, DT_, IT_, 3> o(op);
    TEST_CHECK_EQUAL(o, g);
#ifdef FEAT_HAVE_ZLIB
    auto zb = tester_blocked.serialize(LAFEM::SerialConfig(true,false));
    DenseVectorBlocked<Mem_, DT_, IT_, 5> zlib(zb);
    TEST_CHECK_EQUAL(zlib, tester_blocked);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zp = tester_blocked.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    DenseVectorBlocked<Mem_, DT_, IT_, 5> zfp(zp);
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
DenseVectorBlockedSerializeTest<Mem::Main, float, unsigned int> cpu_dense_vector_blocked_serialize_test_float_uint;
DenseVectorBlockedSerializeTest<Mem::Main, double, unsigned int> cpu_dense_vector_blocked_serialize_test_double_uint;
DenseVectorBlockedSerializeTest<Mem::Main, float, unsigned long> cpu_dense_vector_blocked_serialize_test_float_ulong;
DenseVectorBlockedSerializeTest<Mem::Main, double, unsigned long> cpu_dense_vector_blocked_serialize_test_double_ulong;
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedSerializeTest<Mem::CUDA, float, unsigned int> cuda_dense_vector_blocked_serialize_test_float_uint;
DenseVectorBlockedSerializeTest<Mem::CUDA, double, unsigned int> cuda_dense_vector_blocked_serialize_test_double_uint;
DenseVectorBlockedSerializeTest<Mem::CUDA, float, unsigned long> cuda_dense_vector_blocked_serialize_test_float_ulong;
DenseVectorBlockedSerializeTest<Mem::CUDA, double, unsigned long> cuda_dense_vector_blocked_serialize_test_double_ulong;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedAxpyTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedAxpyTest()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedAxpyTest")
  {
  }

  virtual ~DenseVectorBlockedAxpyTest()
  {
  }

  virtual void run() const override
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> b_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> ref(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j] = DT_(i % 100) * DT_(1.234) + DT_(j);
        a_local(i, tv1);
        Tiny::Vector<DT_, BS_> tv2;
        for (Index j(0) ; j < BS_ ; ++j)
          tv2.v[j] = DT_(2 - DT_(i % (42 + j)));
        b_local(i, tv2);
        ref(i, s * a_local(i) + b_local(i));
      }
      DenseVectorBlocked<Mem_, DT_, IT_, BS_> a(size);
      a.copy(a_local);
      DenseVectorBlocked<Mem_, DT_, IT_, BS_> b(size);
      b.copy(b_local);

      DenseVectorBlocked<Mem_, DT_, IT_, BS_> c(size);
      c.axpy(a, b, s);
      result_local.copy(c);
      for (Index i(0) ; i < size ; ++i)
        for (Index j(0) ; j < BS_ ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i).v[j], ref(i).v[j], DT_(1e-2));

      a.axpy(a, b, s);
      result_local.copy(a);
      for (Index i(0) ; i < size ; ++i)
        for (Index j(0) ; j < BS_ ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i).v[j], ref(i).v[j], DT_(1e-2));

      a.copy(a_local);
      b.axpy(a, b, s);
      result_local.copy(b);
      for (Index i(0) ; i < size ; ++i)
        for (Index j(0) ; j < BS_ ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i).v[j], ref(i).v[j], DT_(1e-2));
    }
  }
};
DenseVectorBlockedAxpyTest<Mem::Main, float, unsigned int, 2> dv_axpy_test_float_uint;
DenseVectorBlockedAxpyTest<Mem::Main, double, unsigned int, 2> dv_axpy_test_double_uint;
DenseVectorBlockedAxpyTest<Mem::Main, float, unsigned long, 3> dv_axpy_test_float_ulong;
DenseVectorBlockedAxpyTest<Mem::Main, double, unsigned long, 3> dv_axpy_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedAxpyTest<Mem::Main, __float128, unsigned int, 2> dv_axpy_test_float128_uint;
DenseVectorBlockedAxpyTest<Mem::Main, __float128, unsigned long, 3> dv_axpy_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedAxpyTest<Mem::CUDA, float, unsigned int, 2> cuda_dv_axpy_test_float_uint;
DenseVectorBlockedAxpyTest<Mem::CUDA, double, unsigned int, 2> cuda_dv_axpy_test_double_uint;
DenseVectorBlockedAxpyTest<Mem::CUDA, float, unsigned long, 3> cuda_dv_axpy_test_float_ulong;
DenseVectorBlockedAxpyTest<Mem::CUDA, double, unsigned long, 3> cuda_dv_axpy_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedDotTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedDotTest()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedDotTest")
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
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> b_local(size);
      const DT_ den(DT_(1) / DT_(size * BS_));
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j] = DT_(i * BS_ + j + 1) * den;
        a_local(i, tv1);
        Tiny::Vector<DT_, BS_> tv2;
        for (Index j(0) ; j < BS_ ; ++j)
          tv2.v[j] = DT_(1) / DT_(i * BS_ + j + 1);
        b_local(i, tv2);
      }

      DenseVectorBlocked<Mem_, DT_, IT_, BS_> a;
      a.convert(a_local);
      DenseVectorBlocked<Mem_, DT_, IT_, BS_> b;
      b.convert(b_local);

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
DenseVectorBlockedDotTest<Mem::Main, float, unsigned int, 3> dv_dot_product_test_float_uint;
DenseVectorBlockedDotTest<Mem::Main, double, unsigned int, 3> dv_dot_product_test_double_uint;
DenseVectorBlockedDotTest<Mem::Main, float, unsigned long, 2> dv_dot_product_test_float_ulong;
DenseVectorBlockedDotTest<Mem::Main, double, unsigned long, 2> dv_dot_product_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedDotTest<Mem::Main, __float128, unsigned int, 3> dv_dot_product_test_float128_uint;
DenseVectorBlockedDotTest<Mem::Main, __float128, unsigned long, 2> dv_dot_product_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedDotTest<Mem::CUDA, float, unsigned int, 3> cuda_dv_dot_product_test_float_uint;
DenseVectorBlockedDotTest<Mem::CUDA, double, unsigned int, 3> cuda_dv_dot_product_test_double_uint;
DenseVectorBlockedDotTest<Mem::CUDA, float, unsigned long, 2> cuda_dv_dot_product_test_float_ulong;
DenseVectorBlockedDotTest<Mem::CUDA, double, unsigned long, 2> cuda_dv_dot_product_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedTripleDotTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedTripleDotTest()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedTripleDotTest")
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
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> b_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> c_local(size);
      const DT_ den(DT_(1) / DT_(size * BS_));
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j] = DT_(i * BS_ + j + 1) * den;
        a_local(i, tv1);
        Tiny::Vector<DT_, BS_> tv2;
        for (Index j(0) ; j < BS_ ; ++j)
          tv2.v[j] = DT_(1) / DT_(i * BS_ + j + 1);
        b_local(i, tv2);
        Tiny::Vector<DT_, BS_> tv3;
        for (Index j(0) ; j < BS_ ; ++j)
          tv3.v[j] = DT_(3) / DT_(i * BS_ + j + 1);
        c_local(i, tv3);
      }

      DenseVectorBlocked<Mem_, DT_, IT_, BS_> a;
      a.convert(a_local);
      DenseVectorBlocked<Mem_, DT_, IT_, BS_> b;
      b.convert(b_local);
      DenseVectorBlocked<Mem_, DT_, IT_, BS_> c;
      c.convert(c_local);

      DenseVector<Mem_, DT_, IT_> ref_a;
      ref_a.convert(a);
      DenseVector<Mem_, DT_, IT_> ref_b;
      ref_b.convert(b);
      DenseVector<Mem_, DT_, IT_> ref_c;
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
DenseVectorBlockedTripleDotTest<Mem::Main, float, unsigned int, 3> dv_triple_dot_product_test_float_uint;
DenseVectorBlockedTripleDotTest<Mem::Main, double, unsigned int, 3> dv_triple_dot_product_test_double_uint;
DenseVectorBlockedTripleDotTest<Mem::Main, float, unsigned long, 2> dv_triple_dot_product_test_float_ulong;
DenseVectorBlockedTripleDotTest<Mem::Main, double, unsigned long, 2> dv_triple_dot_product_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedTripleDotTest<Mem::Main, __float128, unsigned int, 3> dv_triple_dot_product_test_float128_uint;
DenseVectorBlockedTripleDotTest<Mem::Main, __float128, unsigned long, 2> dv_triple_dot_product_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedTripleDotTest<Mem::CUDA, float, unsigned int, 3> cuda_dv_triple_dot_product_test_float_uint;
DenseVectorBlockedTripleDotTest<Mem::CUDA, double, unsigned int, 3> cuda_dv_triple_dot_product_test_double_uint;
DenseVectorBlockedTripleDotTest<Mem::CUDA, float, unsigned long, 2> cuda_dv_triple_dot_product_test_float_ulong;
DenseVectorBlockedTripleDotTest<Mem::CUDA, double, unsigned long, 2> cuda_dv_triple_dot_product_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedComponentProductTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedComponentProductTest()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedComponentProductTest")
  {
  }

  virtual ~DenseVectorBlockedComponentProductTest()
  {
  }

  void run1() const
  {
    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> b_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> ref(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> ref2(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> result_local(size);

      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j] = DT_(DT_(j) + DT_(i) * DT_(1.234));
        a_local(i, tv1);
        Tiny::Vector<DT_, BS_> tv2;
        for (Index j(0) ; j < BS_ ; ++j)
          tv2.v[j] = DT_(size*2*BS_ - i + 2*j);
        b_local(i, tv2);
        ref(i, component_product(a_local(i), b_local(i)) );
        ref2(i, component_product(a_local(i), a_local(i)) );
      }

      DenseVectorBlocked<Mem_, DT_, IT_, BS_> a(size);
      a.copy(a_local);
      DenseVectorBlocked<Mem_, DT_, IT_, BS_> b(size);
      b.copy(b_local);
      DenseVectorBlocked<Mem_, DT_, IT_, BS_> c(size);

      c.component_product(a, b);
      result_local.copy(c);
      TEST_CHECK_EQUAL(result_local, ref);

      a.component_product(a, b);
      result_local.copy(a);
      TEST_CHECK_EQUAL(result_local, ref);

      a.copy(a_local);
      b.component_product(a, b);
      result_local.copy(b);
      TEST_CHECK_EQUAL(result_local, ref);

      b.copy(b_local);
      a.component_product(a, a);
      result_local.copy(a);
      TEST_CHECK_EQUAL(result_local, ref2);
    }
  }

  virtual void run() const override
  {
    run1();
  }
};
DenseVectorBlockedComponentProductTest<Mem::Main, float, unsigned int, 3> dv_component_product_test_float_uint;
DenseVectorBlockedComponentProductTest<Mem::Main, double, unsigned int, 3> dv_component_product_test_double_uint;
DenseVectorBlockedComponentProductTest<Mem::Main, float, unsigned long, 2> dv_component_product_test_float_ulong;
DenseVectorBlockedComponentProductTest<Mem::Main, double, unsigned long, 2> dv_component_product_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedComponentProductTest<Mem::Main, __float128, unsigned int, 3> dv_component_product_test_float128_uint;
DenseVectorBlockedComponentProductTest<Mem::Main, __float128, unsigned long, 2> dv_component_product_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedComponentProductTest<Mem::CUDA, float, unsigned int, 3> cuda_dv_component_product_test_float_uint;
DenseVectorBlockedComponentProductTest<Mem::CUDA, double, unsigned int, 3> cuda_dv_component_product_test_double_uint;
DenseVectorBlockedComponentProductTest<Mem::CUDA, float, unsigned long, 2> cuda_dv_component_product_test_float_ulong;
DenseVectorBlockedComponentProductTest<Mem::CUDA, double, unsigned long, 2> cuda_dv_component_product_test_double_ulong;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedScaleTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedScaleTest()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedScaleTest")
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
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> ref(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j] = DT_(DT_(i) * DT_(BS_) + DT_(j) * DT_(1.234));
        a_local(i, tv1);
        ref(i, a_local(i) * s);
      }

      DenseVectorBlocked<Mem_, DT_, IT_, BS_> a(size);
      a.copy(a_local);
      DenseVectorBlocked<Mem_, DT_, IT_, BS_> b(size);

      b.scale(a, s);
      result_local.copy(b);
      TEST_CHECK_EQUAL(result_local, ref);

      a.scale(a, s);
      result_local.copy(a);
      TEST_CHECK_EQUAL(result_local, ref);
    }
  }
};
DenseVectorBlockedScaleTest<Mem::Main, float, unsigned int, 3> dv_scale_test_float_uint;
DenseVectorBlockedScaleTest<Mem::Main, double, unsigned int, 3> dv_scale_test_double_uint;
DenseVectorBlockedScaleTest<Mem::Main, float, unsigned long, 2> dv_scale_test_float_ulong;
DenseVectorBlockedScaleTest<Mem::Main, double, unsigned long, 2> dv_scale_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedScaleTest<Mem::Main, __float128, unsigned int, 3> dv_scale_test_float128_uint;
DenseVectorBlockedScaleTest<Mem::Main, __float128, unsigned long, 2> dv_scale_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedScaleTest<Mem::CUDA, float, unsigned int, 3> cuda_dv_scale_test_float_uint;
DenseVectorBlockedScaleTest<Mem::CUDA, double, unsigned int, 3> cuda_dv_scale_test_double_uint;
DenseVectorBlockedScaleTest<Mem::CUDA, float, unsigned long, 2> cuda_dv_scale_test_float_ulong;
DenseVectorBlockedScaleTest<Mem::CUDA, double, unsigned long, 2> cuda_dv_scale_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedNorm2Test
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedNorm2Test()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedNorm2Test")
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
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        // a[i] = 1/sqrt(2^i) = (1/2)^(i/2)
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j] = Math::pow(DT_(0.5), DT_(0.5) * DT_(i * BS_ + j));
        a_local(i, tv1);
      }

      // ||a||_2 = sqrt(2 - 2^{1-n})
      const DT_ ref(Math::sqrt(DT_(2) - Math::pow(DT_(0.5), DT_(size*BS_-1))));

      DenseVectorBlocked<Mem_, DT_, IT_, BS_> a;
      a.convert(a_local);
      DT_ c = a.norm2();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);

      c = a.norm2sqr();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref*ref, eps);
    }
  }
};
DenseVectorBlockedNorm2Test<Mem::Main, float, unsigned int, 2> dv_norm2_test_float_uint;
DenseVectorBlockedNorm2Test<Mem::Main, double, unsigned int, 2> dv_norm2_test_double_uint;
DenseVectorBlockedNorm2Test<Mem::Main, float, unsigned long, 3> dv_norm2_test_float_ulong;
DenseVectorBlockedNorm2Test<Mem::Main, double, unsigned long, 3> dv_norm2_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedNorm2Test<Mem::Main, __float128, unsigned int, 2> dv_norm2_test_float128_uint;
DenseVectorBlockedNorm2Test<Mem::Main, __float128, unsigned long, 3> dv_norm2_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedNorm2Test<Mem::CUDA, float, unsigned int, 2> cuda_dv_norm2_test_float_uint;
DenseVectorBlockedNorm2Test<Mem::CUDA, double, unsigned int, 2> cuda_dv_norm2_test_double_uint;
DenseVectorBlockedNorm2Test<Mem::CUDA, float, unsigned long, 3> cuda_dv_norm2_test_float_ulong;
DenseVectorBlockedNorm2Test<Mem::CUDA, double, unsigned long, 3> cuda_dv_norm2_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedMaxAbsElementTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedMaxAbsElementTest()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedMaxAbsElementTest")
  {
  }

  virtual ~DenseVectorBlockedMaxAbsElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e4) ; size*=2)
    {
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j]  = DT_(DT_((i * BS_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        a_local(i, tv1);
      }

      DenseVectorBlocked<Mem_, DT_, IT_, BS_> a;
      a.convert(a_local);
      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ max = a.max_abs_element();

      TEST_CHECK_EQUAL(max, DT_((size*BS_) -1));
    }
  }
};
DenseVectorBlockedMaxAbsElementTest<Mem::Main, float, unsigned int, 2> dv_max_abs_element_test_float_uint;
DenseVectorBlockedMaxAbsElementTest<Mem::Main, double, unsigned int, 2> dv_max_abs_element_test_double_uint;
DenseVectorBlockedMaxAbsElementTest<Mem::Main, float, unsigned long, 3> dv_max_abs_element_test_float_ulong;
DenseVectorBlockedMaxAbsElementTest<Mem::Main, double, unsigned long, 3> dv_max_abs_element_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedMaxAbsElementTest<Mem::Main, __float128, unsigned int, 2> dv_max_abs_element_test_float128_uint;
DenseVectorBlockedMaxAbsElementTest<Mem::Main, __float128, unsigned long, 3> dv_max_abs_element_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedMaxAbsElementTest<Mem::CUDA, float, unsigned int, 2> cuda_dv_max_abs_element_test_float_uint;
DenseVectorBlockedMaxAbsElementTest<Mem::CUDA, double, unsigned int, 2> cuda_dv_max_abs_element_test_double_uint;
DenseVectorBlockedMaxAbsElementTest<Mem::CUDA, float, unsigned long, 3> cuda_dv_max_abs_element_test_float_ulong;
DenseVectorBlockedMaxAbsElementTest<Mem::CUDA, double, unsigned long, 3> cuda_dv_max_abs_element_test_double_ulong;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedMinAbsElementTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedMinAbsElementTest()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedMinAbsElementTest")
  {
  }

  virtual ~DenseVectorBlockedMinAbsElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e4) ; size*=2)
    {
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j]  = DT_(DT_((i * BS_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        a_local(i, tv1);
      }

      DenseVectorBlocked<Mem_, DT_, IT_, BS_> a;
      a.convert(a_local);
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
DenseVectorBlockedMinAbsElementTest<Mem::Main, float, unsigned int, 2> dv_min_abs_element_test_float_uint;
DenseVectorBlockedMinAbsElementTest<Mem::Main, double, unsigned int, 2> dv_min_abs_element_test_double_uint;
DenseVectorBlockedMinAbsElementTest<Mem::Main, float, unsigned long, 3> dv_min_abs_element_test_float_ulong;
DenseVectorBlockedMinAbsElementTest<Mem::Main, double, unsigned long, 3> dv_min_abs_element_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
DenseVectorBlockedMinAbsElementTest<Mem::Main, __float128, unsigned int, 2> dv_min_abs_element_test_float128_uint;
DenseVectorBlockedMinAbsElementTest<Mem::Main, __float128, unsigned long, 3> dv_min_abs_element_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedMinAbsElementTest<Mem::CUDA, float, unsigned int, 2> cuda_dv_min_abs_element_test_float_uint;
DenseVectorBlockedMinAbsElementTest<Mem::CUDA, double, unsigned int, 2> cuda_dv_min_abs_element_test_double_uint;
DenseVectorBlockedMinAbsElementTest<Mem::CUDA, float, unsigned long, 3> cuda_dv_min_abs_element_test_float_ulong;
DenseVectorBlockedMinAbsElementTest<Mem::CUDA, double, unsigned long, 3> cuda_dv_min_abs_element_test_double_ulong;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedPermuteTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  DenseVectorBlockedPermuteTest()
    : FullTaggedTest<Mem_, DT_, IT_>("DenseVectorBlockedPermuteTest")
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
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j]  = DT_(DT_((i * BS_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        a_local(i, tv1);
      }

      DenseVectorBlocked<Mem_, DT_, IT_, BS_> a;
      a.convert(a_local);
      DenseVectorBlocked<Mem_, DT_, IT_, BS_> ref(a.size());
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
DenseVectorBlockedPermuteTest<Mem::Main, float, unsigned int, 2> dv_permute_test_float_uint;
DenseVectorBlockedPermuteTest<Mem::Main, double, unsigned int, 2> dv_permute_test_double_uint;
DenseVectorBlockedPermuteTest<Mem::Main, float, unsigned long, 3> dv_permute_test_float_ulong;
DenseVectorBlockedPermuteTest<Mem::Main, double, unsigned long, 3> dv_permute_test_double_ulong;
#ifdef FEAT_HAVE_CUDA
DenseVectorBlockedPermuteTest<Mem::CUDA, float, unsigned int, 2> cuda_dv_permute_test_float_uint;
DenseVectorBlockedPermuteTest<Mem::CUDA, double, unsigned int, 2> cuda_dv_permute_test_double_uint;
DenseVectorBlockedPermuteTest<Mem::CUDA, float, unsigned long, 3> cuda_dv_permute_test_float_ulong;
DenseVectorBlockedPermuteTest<Mem::CUDA, double, unsigned long, 3> cuda_dv_permute_test_double_ulong;
#endif
