#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/axpy.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class DVAxpyTest
  : public TaggedTest<Arch_, DT_>
{

public:

  DVAxpyTest()
    : TaggedTest<Arch_, DT_>("dv_axpy_test")
  {
  }

  virtual void run() const
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < 1e5 ; size*=2)
    {
      DenseVector<Arch_, DT_> a(size);
      DenseVector<Arch_, DT_> b(size);
      DenseVector<Arch_, DT_> c(size);
      DenseVector<Archs::CPU, DT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i % 100 * DT_(1.234)));
        b(i, DT_(2 - i % 100));
        ref(i, s * a(i) + b(i));
      }

      Axpy<Arch_, BType_>::value(c, s, a, b);
      DenseVector<Archs::CPU, DT_> c2(c);
      for (Index i(0) ; i < c2.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(c2(i), ref(i), 1e-2);

      Axpy<Arch_, BType_>::value(a, s, a, b);
      DenseVector<Archs::CPU, DT_> a2(a);
      for (Index i(0) ; i < a2.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(a2(i), ref(i), 1e-2);
    }
  }
};
DVAxpyTest<Archs::CPU, Archs::Generic, float> dv_axpy_test_float;
DVAxpyTest<Archs::CPU, Archs::Generic, double> dv_axpy_test_double;
#ifdef FEAST_BACKENDS_CUDA
DVAxpyTest<Archs::GPU, Archs::CUDA, float> gpu_dv_axpy_test_float;
DVAxpyTest<Archs::GPU, Archs::CUDA, double> gpu_dv_axpy_test_double;
#endif

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class DVAxpyVTest
  : public TaggedTest<Arch_, DT_>
{

public:

  DVAxpyVTest()
    : TaggedTest<Arch_, DT_>("dv_axpy_v_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e5 ; size*=2)
    {
      DenseVector<Arch_, DT_> a(size);
      DenseVector<Arch_, DT_> b(size);
      DenseVector<Arch_, DT_> c(size);
      DenseVector<Arch_, DT_> d(size);
      DenseVector<Archs::CPU, DT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i % 100 * DT_(1.234)));
        b(i, DT_(2 - i % 100));
        c(i, DT_(1 - i % 100));
        ref(i, c(i) * a(i) + b(i));
      }

      Axpy<Arch_, BType_>::value(d, c, a, b);
      DenseVector<Archs::CPU, DT_> d2(d);
      for (Index i(0) ; i < d2.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(d2(i), ref(i), 1e-2);
      Axpy<Arch_, BType_>::value(a, c, a, b);
      DenseVector<Archs::CPU, DT_> a2(a);
      for (Index i(0) ; i < a2.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(a2(i), ref(i), 1e-2);
    }
  }
};
DVAxpyVTest<Archs::CPU, Archs::Generic, float> dv_axpy_v_test_float;
DVAxpyVTest<Archs::CPU, Archs::Generic, double> dv_axpy_v_test_double;
#ifdef FEAST_BACKENDS_CUDA
DVAxpyVTest<Archs::GPU, Archs::CUDA, float> gpu_dv_axpy_v_test_float;
DVAxpyVTest<Archs::GPU, Archs::CUDA, double> gpu_dv_axpy_v_test_double;
#endif
