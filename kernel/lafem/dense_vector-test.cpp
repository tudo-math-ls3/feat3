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
  typename DT_>
class DenseVectorTest
  : public TaggedTest<Mem_, DT_>
{
public:
  DenseVectorTest()
    : TaggedTest<Mem_, DT_>("DenseVectorTest")
  {
  }

  virtual void run() const
  {
    DenseVector<Mem_, DT_> zero1;
    DenseVector<Mem::Main, DT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    DenseVector<Mem_, DT_> a(10, DT_(7));
    TEST_CHECK_EQUAL(a.bytes_allocated(), 10 * sizeof(DT_) + sizeof(Index));
    DenseVector<Mem_, DT_> b(10, DT_(5));
    b(7, DT_(42));
    DenseVector<Mem_, DT_> c(b.clone());
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

    DenseVector<Mem_, DT_> g(b.size(), b.elements());
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


    DenseVector<Mem_, DT_> k(123);
    for (Index i(0) ; i < k.size() ; ++i)
      k(i, DT_(i) / DT_(12));

    std::stringstream ts;
    k.write_out(FileMode::fm_exp, ts);
    DenseVector<Mem_, DT_> l(FileMode::fm_exp, ts);
    for (Index i(0) ; i < k.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(l(i), k(i), 1e-4);

    BinaryStream bs;
    k.write_out(FileMode::fm_dv, bs);
    bs.seekg(0);
    DenseVector<Mem_, DT_> m(FileMode::fm_dv, bs);
    for (Index i(0) ; i < k.size() ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(m(i), k(i), 1e-5);
  }
};
DenseVectorTest<Mem::Main, float> cpu_dense_vector_test_float;
DenseVectorTest<Mem::Main, double> cpu_dense_vector_test_double;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorTest<Mem::Main, __float128> cpu_dense_vector_test_float128;
#endif
//DenseVectorTest<Mem::Main, Index> cpu_dense_vector_test_index;
#ifdef FEAST_BACKENDS_CUDA
DenseVectorTest<Mem::CUDA, float> cuda_dense_vector_test_float;
DenseVectorTest<Mem::CUDA, double> cuda_dense_vector_test_double;
//DenseVectorTest<Mem::CUDA, Index> cuda_dense_vector_test_index;
#endif

template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class DenseVectorAxpyTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  DenseVectorAxpyTest()
    : TaggedTest<Mem_, DT_, Algo_>("DenseVectorAxpyTest")
  {
  }

  virtual void run() const
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i % 100 * DT_(1.234)));
        b_local(i, DT_(2 - DT_(i % 42)));
        ref(i, s * a_local(i) + b_local(i));
      }
      DenseVector<Mem_, DT_> a(size);
      a.copy(a_local);
      DenseVector<Mem_, DT_> b(size);
      b.copy(b_local);

      DenseVector<Mem_, DT_> c(size);
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
DenseVectorAxpyTest<Mem::Main, Algo::Generic, float> dv_axpy_test_float;
DenseVectorAxpyTest<Mem::Main, Algo::Generic, double> dv_axpy_test_double;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorAxpyTest<Mem::Main, Algo::Generic, __float128> dv_axpy_test_float128;
#endif
#ifdef FEAST_BACKENDS_MKL
DenseVectorAxpyTest<Mem::Main, Algo::MKL, float> mkl_dv_axpy_test_float;
DenseVectorAxpyTest<Mem::Main, Algo::MKL, double> mkl_dv_axpy_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorAxpyTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_axpy_test_float;
DenseVectorAxpyTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_axpy_test_double;
#endif

template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class DenseVectorDotTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  DenseVectorDotTest()
    : TaggedTest<Mem_, DT_, Algo_>("DenseVectorDotTest")
  {
  }

  virtual void run() const
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      const DT_ den(DT_(1) / DT_(size));
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i+1) * den);    // a[i] = (i+1) / n
        b_local(i, DT_(1) / DT_(i+1)); // b[i] = 1 / (i+1)
      }

      DenseVector<Mem_, DT_> a;
      a.convert(a_local);
      DenseVector<Mem_, DT_> b;
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
DenseVectorDotTest<Mem::Main, Algo::Generic, float> dv_dot_product_test_float;
DenseVectorDotTest<Mem::Main, Algo::Generic, double> dv_dot_product_test_double;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorDotTest<Mem::Main, Algo::Generic, __float128> dv_dot_product_test_float128;
#endif
#ifdef FEAST_BACKENDS_MKL
DenseVectorDotTest<Mem::Main, Algo::MKL, float> mkl_dv_dot_product_test_float;
DenseVectorDotTest<Mem::Main, Algo::MKL, double> mkl_dv_dot_product_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorDotTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_dot_product_test_float;
DenseVectorDotTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_dot_product_test_double;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class DenseVectorComponentProductTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  DenseVectorComponentProductTest()
    : TaggedTest<Mem_, DT_, Algo_>("DenseVectorComponentProductTest")
  {
  }

  void run1() const
  {
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> ref2(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i * DT_(1.234)));
        b_local(i, DT_(size*2 - i));
        ref(i, a_local(i) * b_local(i));
        ref2(i, a_local(i) * a_local(i));
      }

      DenseVector<Mem_, DT_> a(size);
      a.copy(a_local);
      DenseVector<Mem_, DT_> b(size);
      b.copy(b_local);
      DenseVector<Mem_, DT_> c(size);

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
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      DenseVector<Mem::Main, DT_> c_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i % 100 * DT_(1.234)));
        b_local(i, DT_(2 - DT_(i % 42)));
        c_local(i, DT_(1 - DT_(i % 23)));
        ref(i, c_local(i) * a_local(i) + b_local(i));
      }
      DenseVector<Mem_, DT_> a(size);
      a.copy(a_local);
      DenseVector<Mem_, DT_> b(size);
      b.copy(b_local);
      DenseVector<Mem_, DT_> c(size);
      c.copy(c_local);

      DenseVector<Mem_, DT_> d(size);
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
DenseVectorComponentProductTest<Mem::Main, Algo::Generic, float> dv_component_product_test_float;
DenseVectorComponentProductTest<Mem::Main, Algo::Generic, double> dv_component_product_test_double;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorComponentProductTest<Mem::Main, Algo::Generic, __float128> dv_component_product_test_float128;
#endif
#ifdef FEAST_BACKENDS_MKL
DenseVectorComponentProductTest<Mem::Main, Algo::MKL, float> mkl_dv_component_product_test_float;
DenseVectorComponentProductTest<Mem::Main, Algo::MKL, double> mkl_dv_component_product_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorComponentProductTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_component_product_test_float;
DenseVectorComponentProductTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_component_product_test_double;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class DenseVectorScaleTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  DenseVectorScaleTest()
    : TaggedTest<Mem_, DT_, Algo_>("DenseVectorScaleTest")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DT_ s(DT_(4.321));
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i * DT_(1.234)));
        ref(i, a_local(i) * s);
      }

      DenseVector<Mem_, DT_> a(size);
      a.copy(a_local);
      DenseVector<Mem_, DT_> b(size);

      b.template scale<Algo_>(a, s);
      result_local.copy(b);
      TEST_CHECK_EQUAL(result_local, ref);

      a.template scale<Algo_>(a, s);
      result_local.copy(a);
      TEST_CHECK_EQUAL(result_local, ref);
    }
  }
};
DenseVectorScaleTest<Mem::Main, Algo::Generic, float> dv_scale_test_float;
DenseVectorScaleTest<Mem::Main, Algo::Generic, double> dv_scale_test_double;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorScaleTest<Mem::Main, Algo::Generic, __float128> dv_scale_test_float128;
#endif
#ifdef FEAST_BACKENDS_MKL
DenseVectorScaleTest<Mem::Main, Algo::MKL, float> mkl_dv_scale_test_float;
DenseVectorScaleTest<Mem::Main, Algo::MKL, double> mkl_dv_scale_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorScaleTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_scale_test_float;
DenseVectorScaleTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_scale_test_double;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class DenseVectorNorm2Test
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  DenseVectorNorm2Test()
    : TaggedTest<Mem_, DT_, Algo_>("DenseVectorNorm2Test")
  {
  }

  virtual void run() const
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_> a_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        // a[i] = 1/sqrt(2^i) = (1/2)^(i/2)
        a_local(i, Math::pow(DT_(0.5), DT_(0.5) * DT_(i)));
      }

      // ||a||_2 = sqrt(2 - 2^{1-n})
      const DT_ ref(Math::sqrt(DT_(2) - Math::pow(DT_(0.5), DT_(size-1))));

      DenseVector<Mem_, DT_> a;
      a.convert(a_local);
      DT_ c = a.template norm2<Algo_>();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);

      c = a.template norm2sqr<Algo_>();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref*ref, eps);
    }
  }
};
DenseVectorNorm2Test<Mem::Main, Algo::Generic, float> dv_norm2_test_float;
DenseVectorNorm2Test<Mem::Main, Algo::Generic, double> dv_norm2_test_double;
#ifdef FEAST_HAVE_QUADMATH
DenseVectorNorm2Test<Mem::Main, Algo::Generic, __float128> dv_norm2_test_float128;
#endif
#ifdef FEAST_BACKENDS_MKL
DenseVectorNorm2Test<Mem::Main, Algo::MKL, float> mkl_dv_norm2_test_float;
DenseVectorNorm2Test<Mem::Main, Algo::MKL, double> mkl_dv_norm2_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DenseVectorNorm2Test<Mem::CUDA, Algo::CUDA, float> cuda_dv_norm2_test_float;
DenseVectorNorm2Test<Mem::CUDA, Algo::CUDA, double> cuda_dv_norm2_test_double;
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
