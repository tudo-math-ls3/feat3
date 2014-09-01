#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/binary_stream.hpp>

#include <list>
#include <sstream>
#include <cstdio>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the dense vector class.
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
class DenseVectorTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorTest")
  {
  }

  virtual void run() const
  {
    DenseVector<Mem_, DT_, IT_> zero1;
    DenseVector<Mem::Main, DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    DenseVector<Mem_, DT_, IT_> a(10, DT_(7));
    TEST_CHECK_EQUAL(a.bytes_allocated(), 10 * sizeof(DT_) + sizeof(Index));
    DenseVector<Mem_, DT_, IT_> b(10, DT_(5));
    b(7, DT_(42));
    TEST_CHECK_EQUAL(b(7), DT_(42));
    TEST_CHECK_EQUAL(b(3), DT_(5));
    DenseVector<Mem_, DT_, IT_> c(b.clone());
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c(7), b(7));
    TEST_CHECK_EQUAL(c, b);
    c.convert(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c(7), b(7));
    TEST_CHECK_EQUAL(c, b);
    DenseVector<Mem::Main, float, unsigned int> d;
    d.convert(c);
    DenseVector<Mem::Main, float, unsigned int> e;
    e.convert(b);
    TEST_CHECK_EQUAL(e.size(), d.size());
    TEST_CHECK_EQUAL(e(7), d(7));
    TEST_CHECK_EQUAL(e, d);

    b.clone(a);
    TEST_CHECK_NOT_EQUAL((void*)b.elements(), (void*)a.elements());
    c.convert(a);
    TEST_CHECK_EQUAL((void*)c.elements(), (void*)a.elements());
    TEST_CHECK_EQUAL(b, c);
    a(3, DT_(23));
    TEST_CHECK_EQUAL(a, c);
    TEST_CHECK_NOT_EQUAL(a, b);

    DenseVector<Mem_, DT_, IT_> g(b.size(), b.elements());
    TEST_CHECK_EQUAL(g, b);
    TEST_CHECK_EQUAL((void*)g.elements(), (void*)b.elements());

    {
      EDI<Mem_, DT_> t(a.edi(2));
      t = DT_(41);
      TEST_CHECK_NOT_EQUAL(a(2), DT_(41));
      a.edi(1) = DT_(4);
      TEST_CHECK_EQUAL(a(1), DT_(4));
      a.edi(1) += DT_(4);
      TEST_CHECK_EQUAL(a(1), DT_(8));
    }
    TEST_CHECK_EQUAL(a(1), DT_(8));
    TEST_CHECK_EQUAL(a(2), DT_(41));


    DenseVector<Mem_, DT_, IT_> k(123);
    for (Index i(0) ; i < k.size() ; ++i)
      k(i, DT_(i) / DT_(12));

    std::stringstream mts;
    k.write_out(FileMode::fm_mtx, mts);
    DenseVector<Mem_, DT_, IT_> l(FileMode::fm_mtx, mts);
    for (Index i(0) ; i < k.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(l(i), k(i), 1e-4);

    std::stringstream ts;
    k.write_out(FileMode::fm_exp, ts);
    DenseVector<Mem_, DT_, IT_> m(FileMode::fm_exp, ts);
    for (Index i(0) ; i < k.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(m(i), k(i), 1e-4);

    BinaryStream bs;
    k.write_out(FileMode::fm_dv, bs);
    bs.seekg(0);
    DenseVector<Mem_, DT_, IT_> n(FileMode::fm_dv, bs);
    for (Index i(0) ; i < k.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(n(i), k(i), 1e-5);

    auto op = k.serialise();
    DenseVector<Mem_, DT_, IT_> o(op);
    delete[] op.second;
    for (Index i(0) ; i < k.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(o(i), k(i), 1e-5);
  }
};
DenseVectorTest<Mem::Main, NotSet, float, unsigned int> cpu_dense_vector_test_float_uint;
DenseVectorTest<Mem::Main, NotSet, double, unsigned int> cpu_dense_vector_test_double_uint;
DenseVectorTest<Mem::Main, NotSet, float, unsigned long> cpu_dense_vector_test_float_ulong;
DenseVectorTest<Mem::Main, NotSet, double, unsigned long> cpu_dense_vector_test_double_ulong;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorTest<Mem::Main, NotSet, __float128, unsigned int> cpu_dense_vector_test_float128_uint;
DenseVectorTest<Mem::Main, NotSet, __float128, unsigned long> cpu_dense_vector_test_float128_ulong;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorTest<Mem::CUDA, NotSet, float, unsigned int> cuda_dense_vector_test_float_uint;
DenseVectorTest<Mem::CUDA, NotSet, double, unsigned int> cuda_dense_vector_test_double_uint;
DenseVectorTest<Mem::CUDA, NotSet, float, unsigned long> cuda_dense_vector_test_float_ulong;
DenseVectorTest<Mem::CUDA, NotSet, double, unsigned long> cuda_dense_vector_test_double_ulong;
#endif

template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class DenseVectorAxpyTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorAxpyTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorAxpyTest")
  {
  }

  virtual void run() const
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_, IT_> a_local(size);
      DenseVector<Mem::Main, DT_, IT_> b_local(size);
      DenseVector<Mem::Main, DT_, IT_> ref(size);
      DenseVector<Mem::Main, DT_, IT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i % 100 * DT_(1.234)));
        b_local(i, DT_(2 - DT_(i % 42)));
        ref(i, s * a_local(i) + b_local(i));
      }
      DenseVector<Mem_, DT_, IT_> a(size);
      a.copy(a_local);
      DenseVector<Mem_, DT_, IT_> b(size);
      b.copy(b_local);

      DenseVector<Mem_, DT_, IT_> c(size);
      c.template axpy<Algo_>(a, b, s);
      result_local.copy(c);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      a.template axpy<Algo_>(a, b, s);
      result_local.copy(a);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      a.copy(a_local);
      b.template axpy<Algo_>(a, b, s);
      result_local.copy(b);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);
    }
  }
};
DenseVectorAxpyTest<Mem::Main, Algo::Generic, float, unsigned int> dv_axpy_test_float_uint;
DenseVectorAxpyTest<Mem::Main, Algo::Generic, double, unsigned int> dv_axpy_test_double_uint;
DenseVectorAxpyTest<Mem::Main, Algo::Generic, float, unsigned long> dv_axpy_test_float_ulong;
DenseVectorAxpyTest<Mem::Main, Algo::Generic, double, unsigned long> dv_axpy_test_double_ulong;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorAxpyTest<Mem::Main, Algo::Generic, __float128, unsigned int> dv_axpy_test_float128_uint;
DenseVectorAxpyTest<Mem::Main, Algo::Generic, __float128, unsigned long> dv_axpy_test_float128_ulong;
#endif
#ifdef FEAST_BACKENDS_MKL
DenseVectorAxpyTest<Mem::Main, Algo::MKL, float, unsigned int> mkl_dv_axpy_test_float_uint;
DenseVectorAxpyTest<Mem::Main, Algo::MKL, double, unsigned int> mkl_dv_axpy_test_double_uint;
DenseVectorAxpyTest<Mem::Main, Algo::MKL, float, unsigned long> mkl_dv_axpy_test_float_ulong;
DenseVectorAxpyTest<Mem::Main, Algo::MKL, double, unsigned long> mkl_dv_axpy_test_double_ulong;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorAxpyTest<Mem::CUDA, Algo::CUDA, float, unsigned int> cuda_dv_axpy_test_float_uint;
DenseVectorAxpyTest<Mem::CUDA, Algo::CUDA, double, unsigned int> cuda_dv_axpy_test_double_uint;
DenseVectorAxpyTest<Mem::CUDA, Algo::CUDA, float, unsigned long> cuda_dv_axpy_test_float_ulong;
DenseVectorAxpyTest<Mem::CUDA, Algo::CUDA, double, unsigned long> cuda_dv_axpy_test_double_ulong;
#endif

template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class DenseVectorDotTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorDotTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorDotTest")
  {
  }

  virtual void run() const
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_, IT_> a_local(size);
      DenseVector<Mem::Main, DT_, IT_> b_local(size);
      const DT_ den(DT_(1) / DT_(size));
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i+1) * den);    // a[i] = (i+1) / n
        b_local(i, DT_(1) / DT_(i+1)); // b[i] = 1 / (i+1)
      }

      DenseVector<Mem_, DT_, IT_> a;
      a.convert(a_local);
      DenseVector<Mem_, DT_, IT_> b;
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
DenseVectorDotTest<Mem::Main, Algo::Generic, float, unsigned int> dv_dot_product_test_float_uint;
DenseVectorDotTest<Mem::Main, Algo::Generic, double, unsigned int> dv_dot_product_test_double_uint;
DenseVectorDotTest<Mem::Main, Algo::Generic, float, unsigned long> dv_dot_product_test_float_ulong;
DenseVectorDotTest<Mem::Main, Algo::Generic, double, unsigned long> dv_dot_product_test_double_ulong;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorDotTest<Mem::Main, Algo::Generic, __float128, unsigned int> dv_dot_product_test_float128_uint;
DenseVectorDotTest<Mem::Main, Algo::Generic, __float128, unsigned long> dv_dot_product_test_float128_ulong;
#endif
#ifdef FEAST_BACKENDS_MKL
DenseVectorDotTest<Mem::Main, Algo::MKL, float, unsigned int> mkl_dv_dot_product_test_float_uint;
DenseVectorDotTest<Mem::Main, Algo::MKL, double, unsigned int> mkl_dv_dot_product_test_double_uint;
DenseVectorDotTest<Mem::Main, Algo::MKL, float, unsigned long> mkl_dv_dot_product_test_float_ulong;
DenseVectorDotTest<Mem::Main, Algo::MKL, double, unsigned long> mkl_dv_dot_product_test_double_ulong;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorDotTest<Mem::CUDA, Algo::CUDA, float, unsigned int> cuda_dv_dot_product_test_float_uint;
DenseVectorDotTest<Mem::CUDA, Algo::CUDA, double, unsigned int> cuda_dv_dot_product_test_double_uint;
DenseVectorDotTest<Mem::CUDA, Algo::CUDA, float, unsigned long> cuda_dv_dot_product_test_float_ulong;
DenseVectorDotTest<Mem::CUDA, Algo::CUDA, double, unsigned long> cuda_dv_dot_product_test_double_ulong;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class DenseVectorComponentProductTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorComponentProductTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorComponentProductTest")
  {
  }

  void run1() const
  {
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_, IT_> a_local(size);
      DenseVector<Mem::Main, DT_, IT_> b_local(size);
      DenseVector<Mem::Main, DT_, IT_> ref(size);
      DenseVector<Mem::Main, DT_, IT_> ref2(size);
      DenseVector<Mem::Main, DT_, IT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i * DT_(1.234)));
        b_local(i, DT_(size*2 - i));
        ref(i, a_local(i) * b_local(i));
        ref2(i, a_local(i) * a_local(i));
      }

      DenseVector<Mem_, DT_, IT_> a(size);
      a.copy(a_local);
      DenseVector<Mem_, DT_, IT_> b(size);
      b.copy(b_local);
      DenseVector<Mem_, DT_, IT_> c(size);

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

  void run2() const
  {
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_, IT_> a_local(size);
      DenseVector<Mem::Main, DT_, IT_> b_local(size);
      DenseVector<Mem::Main, DT_, IT_> c_local(size);
      DenseVector<Mem::Main, DT_, IT_> ref(size);
      DenseVector<Mem::Main, DT_, IT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i % 100 * DT_(1.234)));
        b_local(i, DT_(2 - DT_(i % 42)));
        c_local(i, DT_(1 - DT_(i % 23)));
        ref(i, c_local(i) * a_local(i) + b_local(i));
      }
      DenseVector<Mem_, DT_, IT_> a(size);
      a.copy(a_local);
      DenseVector<Mem_, DT_, IT_> b(size);
      b.copy(b_local);
      DenseVector<Mem_, DT_, IT_> c(size);
      c.copy(c_local);

      DenseVector<Mem_, DT_, IT_> d(size);
      d.template component_product<Algo_>(c, a, b);
      result_local.copy(d);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      a.template component_product<Algo_>(c, a, b);
      result_local.copy(a);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      a.copy(a_local);
      b.template component_product<Algo_>(c, a, b);
      result_local.copy(b);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      b.copy(b_local);
      c.template component_product<Algo_>(c, a, b);
      result_local.copy(c);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);
    }
  }

  virtual void run() const
  {
    run1();
    run2();
  }
};
DenseVectorComponentProductTest<Mem::Main, Algo::Generic, float, unsigned int> dv_component_product_test_float_uint;
DenseVectorComponentProductTest<Mem::Main, Algo::Generic, double, unsigned int> dv_component_product_test_double_uint;
DenseVectorComponentProductTest<Mem::Main, Algo::Generic, float, unsigned long> dv_component_product_test_float_ulong;
DenseVectorComponentProductTest<Mem::Main, Algo::Generic, double, unsigned long> dv_component_product_test_double_ulong;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorComponentProductTest<Mem::Main, Algo::Generic, __float128, unsigned int> dv_component_product_test_float128_uint;
DenseVectorComponentProductTest<Mem::Main, Algo::Generic, __float128, unsigned long> dv_component_product_test_float128_ulong;
#endif
#ifdef FEAST_BACKENDS_MKL
DenseVectorComponentProductTest<Mem::Main, Algo::MKL, float, unsigned int> mkl_dv_component_product_test_float_uint;
DenseVectorComponentProductTest<Mem::Main, Algo::MKL, double, unsigned int> mkl_dv_component_product_test_double_uint;
DenseVectorComponentProductTest<Mem::Main, Algo::MKL, float, unsigned long> mkl_dv_component_product_test_float_ulong;
DenseVectorComponentProductTest<Mem::Main, Algo::MKL, double, unsigned long> mkl_dv_component_product_test_double_ulong;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorComponentProductTest<Mem::CUDA, Algo::CUDA, float, unsigned int> cuda_dv_component_product_test_float_uint;
DenseVectorComponentProductTest<Mem::CUDA, Algo::CUDA, double, unsigned int> cuda_dv_component_product_test_double_uint;
DenseVectorComponentProductTest<Mem::CUDA, Algo::CUDA, float, unsigned long> cuda_dv_component_product_test_float_ulong;
DenseVectorComponentProductTest<Mem::CUDA, Algo::CUDA, double, unsigned long> cuda_dv_component_product_test_double_ulong;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class DenseVectorScaleTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorScaleTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorScaleTest")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DT_ s(DT_(4.321));
      DenseVector<Mem::Main, DT_, IT_> a_local(size);
      DenseVector<Mem::Main, DT_, IT_> ref(size);
      DenseVector<Mem::Main, DT_, IT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i * DT_(1.234)));
        ref(i, a_local(i) * s);
      }

      DenseVector<Mem_, DT_, IT_> a(size);
      a.copy(a_local);
      DenseVector<Mem_, DT_, IT_> b(size);

      b.template scale<Algo_>(a, s);
      result_local.copy(b);
      TEST_CHECK_EQUAL(result_local, ref);

      a.template scale<Algo_>(a, s);
      result_local.copy(a);
      TEST_CHECK_EQUAL(result_local, ref);
    }
  }
};
DenseVectorScaleTest<Mem::Main, Algo::Generic, float, unsigned int> dv_scale_test_float_uint;
DenseVectorScaleTest<Mem::Main, Algo::Generic, double, unsigned int> dv_scale_test_double_uint;
DenseVectorScaleTest<Mem::Main, Algo::Generic, float, unsigned long> dv_scale_test_float_ulong;
DenseVectorScaleTest<Mem::Main, Algo::Generic, double, unsigned long> dv_scale_test_double_ulong;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorScaleTest<Mem::Main, Algo::Generic, __float128, unsigned int> dv_scale_test_float128_uint;
DenseVectorScaleTest<Mem::Main, Algo::Generic, __float128, unsigned long> dv_scale_test_float128_ulong;
#endif
#ifdef FEAST_BACKENDS_MKL
DenseVectorScaleTest<Mem::Main, Algo::MKL, float, unsigned int> mkl_dv_scale_test_float_uint;
DenseVectorScaleTest<Mem::Main, Algo::MKL, double, unsigned int> mkl_dv_scale_test_double_uint;
DenseVectorScaleTest<Mem::Main, Algo::MKL, float, unsigned long> mkl_dv_scale_test_float_ulong;
DenseVectorScaleTest<Mem::Main, Algo::MKL, double, unsigned long> mkl_dv_scale_test_double_ulong;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorScaleTest<Mem::CUDA, Algo::CUDA, float, unsigned int> cuda_dv_scale_test_float_uint;
DenseVectorScaleTest<Mem::CUDA, Algo::CUDA, double, unsigned int> cuda_dv_scale_test_double_uint;
DenseVectorScaleTest<Mem::CUDA, Algo::CUDA, float, unsigned long> cuda_dv_scale_test_float_ulong;
DenseVectorScaleTest<Mem::CUDA, Algo::CUDA, double, unsigned long> cuda_dv_scale_test_double_ulong;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class DenseVectorNorm2Test
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  DenseVectorNorm2Test()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("DenseVectorNorm2Test")
  {
  }

  virtual void run() const
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_, IT_> a_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        // a[i] = 1/sqrt(2^i) = (1/2)^(i/2)
        a_local(i, Math::pow(DT_(0.5), DT_(0.5) * DT_(i)));
      }

      // ||a||_2 = sqrt(2 - 2^{1-n})
      const DT_ ref(Math::sqrt(DT_(2) - Math::pow(DT_(0.5), DT_(size-1))));

      DenseVector<Mem_, DT_, IT_> a;
      a.convert(a_local);
      DT_ c = a.template norm2<Algo_>();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);

      c = a.template norm2sqr<Algo_>();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref*ref, eps);
    }
  }
};
DenseVectorNorm2Test<Mem::Main, Algo::Generic, float, unsigned int> dv_norm2_test_float_uint;
DenseVectorNorm2Test<Mem::Main, Algo::Generic, double, unsigned int> dv_norm2_test_double_uint;
DenseVectorNorm2Test<Mem::Main, Algo::Generic, float, unsigned long> dv_norm2_test_float_ulong;
DenseVectorNorm2Test<Mem::Main, Algo::Generic, double, unsigned long> dv_norm2_test_double_ulong;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorNorm2Test<Mem::Main, Algo::Generic, __float128, unsigned int> dv_norm2_test_float128_uint;
DenseVectorNorm2Test<Mem::Main, Algo::Generic, __float128, unsigned long> dv_norm2_test_float128_ulong;
#endif
#ifdef FEAST_BACKENDS_MKL
DenseVectorNorm2Test<Mem::Main, Algo::MKL, float, unsigned int> mkl_dv_norm2_test_float_uint;
DenseVectorNorm2Test<Mem::Main, Algo::MKL, double, unsigned int> mkl_dv_norm2_test_double_uint;
DenseVectorNorm2Test<Mem::Main, Algo::MKL, float, unsigned long> mkl_dv_norm2_test_float_ulong;
DenseVectorNorm2Test<Mem::Main, Algo::MKL, double, unsigned long> mkl_dv_norm2_test_double_ulong;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorNorm2Test<Mem::CUDA, Algo::CUDA, float, unsigned int> cuda_dv_norm2_test_float_uint;
DenseVectorNorm2Test<Mem::CUDA, Algo::CUDA, double, unsigned int> cuda_dv_norm2_test_double_uint;
DenseVectorNorm2Test<Mem::CUDA, Algo::CUDA, float, unsigned long> cuda_dv_norm2_test_float_ulong;
DenseVectorNorm2Test<Mem::CUDA, Algo::CUDA, double, unsigned long> cuda_dv_norm2_test_double_ulong;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class DenseVectorComponentInvertTest
  : public FullTaggedTest<Mem_, DT_, Algo_, IT_>
{
public:
  typedef DenseVector<Mem_, DT_, IT_> VectorType;

  DenseVectorComponentInvertTest()
    : FullTaggedTest<Mem_, DT_, Algo_, IT_>("DenseVectorComponentInvertTest")
  {
  }

  virtual void run() const
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    const DT_ alpha(Math::pi<DT_>());

    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      // create a vector
      DenseVector<Mem::Main, DT_, IT_> tvec(size);
      for (Index i(0); i < size; ++i)
      {
        tvec(i, DT_(7.63) * DT_(i % 3 + 1) - DT_(9.3));
      }
      VectorType vec;
      vec.convert(tvec);

      VectorType vec2(vec.clone());
      vec2.template component_invert<Algo_>(vec2, alpha);
      vec2.template component_product<Algo_>(vec2, vec);
      for (Index i(0); i < size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(vec2(i), alpha, eps);
      }

      VectorType vec3(size);
      vec3.template component_invert<Algo_>(vec);
      for (Index i(0); i < size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(vec3(i), DT_(1.0) / vec(i), eps);
      }
    }
  }
};

DenseVectorComponentInvertTest<Mem::Main, Algo::Generic, float, Index> dv_component_invert_test_float;
DenseVectorComponentInvertTest<Mem::Main, Algo::Generic, double, Index> dv_component_invert_test_double;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorComponentInvertTest<Mem::Main, Algo::Generic, __float128, Index> dv_component_invert_test_float128;
#endif
#ifdef FEAST_BACKENDS_MKL
//DenseVectorComponentInvertTest<Mem::Main, Algo::MKL, float, Index> mkl_dv_component_invert_test_float;
//DenseVectorComponentInvertTest<Mem::Main, Algo::MKL, double, Index> mkl_dv_component_invert_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorComponentInvertTest<Mem::CUDA, Algo::CUDA, float, Index> cuda_dv_component_invert_test_float;
DenseVectorComponentInvertTest<Mem::CUDA, Algo::CUDA, double, Index> cuda_dv_component_invert_test_double;
#endif
