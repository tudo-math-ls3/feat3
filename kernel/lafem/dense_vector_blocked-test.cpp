#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/binary_stream.hpp>

#include <list>
#include <sstream>
#include <cstdio>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

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
  typename Algo_,
  typename DT_,
  typename IT_>
class DenseVectorBlockedTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorBlockedTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorBlockedTest")
  {
  }

  virtual void run() const
  {
    DenseVectorBlocked<Mem_, DT_, IT_, 2> zero1;
    DenseVectorBlocked<Mem::Main, DT_, IT_, 2> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    DenseVectorBlocked<Mem_, DT_, IT_, 2> a(10, DT_(7));
    TEST_CHECK_EQUAL(a, a);
    TEST_CHECK_EQUAL(a.bytes_allocated(), 20 * sizeof(DT_) + sizeof(Index));
    DenseVectorBlocked<Mem_, DT_, IT_, 2> b(10, DT_(5));
    Tiny::Vector<DT_, 2> tv(42);
    b(7, tv);
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
    c = a.shared();
    TEST_CHECK_EQUAL((void*)c.elements(), (void*)a.elements());
    TEST_CHECK_EQUAL(a, c);

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

    DenseVectorBlocked<Mem_, DT_, IT_, 3> g(f.size(), f.raw_elements());
    TEST_CHECK_EQUAL(g, f);
    TEST_CHECK_EQUAL((void*)g.raw_elements(), (void*)f.raw_elements());
  }
};
DenseVectorBlockedTest<Mem::Main, NotSet, float, Index> cpu_dense_vector_blocked_test_float;
DenseVectorBlockedTest<Mem::Main, NotSet, double, Index> cpu_dense_vector_blocked_test_double;
#ifdef FEAST_BACKENDS_CUDA
DenseVectorBlockedTest<Mem::CUDA, NotSet, float, Index> cuda_dense_vector_blocked_test_float;
DenseVectorBlockedTest<Mem::CUDA, NotSet, double, Index> cuda_dense_vector_blocked_test_double;
#endif

template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedAxpyTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorBlockedAxpyTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorBlockedAxpyTest")
  {
  }

  virtual void run() const
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> b_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> ref(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j] = DT_(i % 100 * DT_(1.234 + j));
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
      c.template axpy<Algo_>(a, b, s);
      result_local.copy(c);
      for (Index i(0) ; i < size ; ++i)
        for (Index j(0) ; j < BS_ ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i).v[j], ref(i).v[j], 1e-2);

      a.template axpy<Algo_>(a, b, s);
      result_local.copy(a);
      for (Index i(0) ; i < size ; ++i)
        for (Index j(0) ; j < BS_ ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i).v[j], ref(i).v[j], 1e-2);

      a.copy(a_local);
      b.template axpy<Algo_>(a, b, s);
      result_local.copy(b);
      for (Index i(0) ; i < size ; ++i)
        for (Index j(0) ; j < BS_ ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i).v[j], ref(i).v[j], 1e-2);
    }
  }
};
DenseVectorBlockedAxpyTest<Mem::Main, Algo::Generic, float, Index, 2> dv_axpy_test_float;
DenseVectorBlockedAxpyTest<Mem::Main, Algo::Generic, double, Index, 2> dv_axpy_test_double;
#ifdef FEAST_BACKENDS_MKL
DenseVectorBlockedAxpyTest<Mem::Main, Algo::MKL, float, Index, 2> mkl_dv_axpy_test_float;
DenseVectorBlockedAxpyTest<Mem::Main, Algo::MKL, double, Index, 2> mkl_dv_axpy_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorBlockedAxpyTest<Mem::CUDA, Algo::CUDA, float, Index, 2> cuda_dv_axpy_test_float;
DenseVectorBlockedAxpyTest<Mem::CUDA, Algo::CUDA, double, Index, 2> cuda_dv_axpy_test_double;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedDotTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorBlockedDotTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorBlockedDotTest")
  {
  }

  virtual void run() const
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(1) ; size < 1e3 ; size*=2)
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
      DT_ c  = a.template dot<Algo_>(b);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = b.template dot<Algo_>(a);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = b.template dot<Algo_>(b);
      ref = b.template norm2<Algo_>();
      ref *= ref;
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
    }
  }
};
DenseVectorBlockedDotTest<Mem::Main, Algo::Generic, float, Index, 2> dv_dot_product_test_float;
DenseVectorBlockedDotTest<Mem::Main, Algo::Generic, double, Index, 2> dv_dot_product_test_double;
#ifdef FEAST_BACKENDS_MKL
DenseVectorBlockedDotTest<Mem::Main, Algo::MKL, float, Index, 2> mkl_dv_dot_product_test_float;
DenseVectorBlockedDotTest<Mem::Main, Algo::MKL, double, Index, 2> mkl_dv_dot_product_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorBlockedDotTest<Mem::CUDA, Algo::CUDA, float, Index, 2> cuda_dv_dot_product_test_float;
DenseVectorBlockedDotTest<Mem::CUDA, Algo::CUDA, double, Index, 2> cuda_dv_dot_product_test_double;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedComponentProductTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorBlockedComponentProductTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorBlockedComponentProductTest")
  {
  }

  void run1() const
  {
    for (Index size(1) ; size < 1e3 ; size*=2)
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
          tv1.v[j] = DT_(j + i  *DT_(1.234));
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

      c.template component_product<Algo_>(a, b);
      result_local.copy(c);
      TEST_CHECK_EQUAL(result_local, ref);

      a.template component_product<Algo_>(a, b);
      result_local.copy(a);
      TEST_CHECK_EQUAL(result_local, ref);

      a.copy(a_local);
      b.template component_product<Algo_>(a, b);
      result_local.copy(b);
      TEST_CHECK_EQUAL(result_local, ref);

      b.copy(b_local);
      a.template component_product<Algo_>(a, a);
      result_local.copy(a);
      TEST_CHECK_EQUAL(result_local, ref2);
    }
  }

  virtual void run() const
  {
    run1();
  }
};
DenseVectorBlockedComponentProductTest<Mem::Main, Algo::Generic, float, Index, 2> dv_component_product_test_float;
DenseVectorBlockedComponentProductTest<Mem::Main, Algo::Generic, double, Index, 2> dv_component_product_test_double;
#ifdef FEAST_BACKENDS_MKL
DenseVectorBlockedComponentProductTest<Mem::Main, Algo::MKL, float, Index, 2> mkl_dv_component_product_test_float;
DenseVectorBlockedComponentProductTest<Mem::Main, Algo::MKL, double, Index, 2> mkl_dv_component_product_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorBlockedComponentProductTest<Mem::CUDA, Algo::CUDA, float, Index, 2> cuda_dv_component_product_test_float;
DenseVectorBlockedComponentProductTest<Mem::CUDA, Algo::CUDA, double, Index, 2> cuda_dv_component_product_test_double;
#endif

template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedScaleTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorBlockedScaleTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorBlockedScaleTest")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DT_ s(DT_(4.321));
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> a_local(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> ref(size);
      DenseVectorBlocked<Mem::Main, DT_, IT_, BS_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        Tiny::Vector<DT_, BS_> tv1;
        for (Index j(0) ; j < BS_ ; ++j)
          tv1.v[j] = DT_(i * BS_ + j * DT_(1.234));
        a_local(i, tv1);
        ref(i, a_local(i) * s);
      }

      DenseVectorBlocked<Mem_, DT_, IT_, BS_> a(size);
      a.copy(a_local);
      DenseVectorBlocked<Mem_, DT_, IT_, BS_> b(size);

      b.template scale<Algo_>(a, s);
      result_local.copy(b);
      TEST_CHECK_EQUAL(result_local, ref);

      a.template scale<Algo_>(a, s);
      result_local.copy(a);
      TEST_CHECK_EQUAL(result_local, ref);
    }
  }
};
DenseVectorBlockedScaleTest<Mem::Main, Algo::Generic, float, Index, 2> dv_scale_test_float;
DenseVectorBlockedScaleTest<Mem::Main, Algo::Generic, double, Index, 2> dv_scale_test_double;
#ifdef FEAST_BACKENDS_MKL
DenseVectorBlockedScaleTest<Mem::Main, Algo::MKL, float, Index, 2> mkl_dv_scale_test_float;
DenseVectorBlockedScaleTest<Mem::Main, Algo::MKL, double, Index, 2> mkl_dv_scale_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorBlockedScaleTest<Mem::CUDA, Algo::CUDA, float, Index, 2> cuda_dv_scale_test_float;
DenseVectorBlockedScaleTest<Mem::CUDA, Algo::CUDA, double, Index, 2> cuda_dv_scale_test_double;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_,
  Index BS_>
class DenseVectorBlockedNorm2Test
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorBlockedNorm2Test()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorBlockedNorm2Test")
  {
  }

  virtual void run() const
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(1) ; size < 1e3 ; size*=2)
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
      DT_ c = a.template norm2<Algo_>();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);

      c = a.template norm2sqr<Algo_>();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref*ref, eps);
    }
  }
};
DenseVectorBlockedNorm2Test<Mem::Main, Algo::Generic, float, Index, 2> dv_norm2_test_float;
DenseVectorBlockedNorm2Test<Mem::Main, Algo::Generic, double, Index, 2> dv_norm2_test_double;
#ifdef FEAST_BACKENDS_MKL
DenseVectorBlockedNorm2Test<Mem::Main, Algo::MKL, float, Index, 2> mkl_dv_norm2_test_float;
DenseVectorBlockedNorm2Test<Mem::Main, Algo::MKL, double, Index, 2> mkl_dv_norm2_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorBlockedNorm2Test<Mem::CUDA, Algo::CUDA, float, Index, 2> cuda_dv_norm2_test_float;
DenseVectorBlockedNorm2Test<Mem::CUDA, Algo::CUDA, double, Index, 2> cuda_dv_norm2_test_double;
#endif
