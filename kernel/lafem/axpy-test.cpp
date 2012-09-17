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
      DenseVector<Archs::CPU, DT_> a_local(size);
      DenseVector<Archs::CPU, DT_> b_local(size);
      DenseVector<Arch_, DT_> c(size);
      DenseVector<Archs::CPU, DT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i % 100 * DT_(1.234)));
        b_local(i, DT_(2 - i % 100));
        ref(i, s * a_local(i) + b_local(i));
      }
      DenseVector<Arch_, DT_> a(a_local);
      DenseVector<Arch_, DT_> b(b_local);

      Axpy<Arch_, BType_>::value(c, s, a, b);
      DenseVector<Archs::CPU, DT_> c_local(c);
      for (Index i(0) ; i < c_local.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(c_local(i), ref(i), 1e-2);

      Axpy<Arch_, BType_>::value(a, s, a, b);
      a_local = a;
      for (Index i(0) ; i < a_local.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(a_local(i), ref(i), 1e-2);
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
      DenseVector<Archs::CPU, DT_> a_local(size);
      DenseVector<Archs::CPU, DT_> b_local(size);
      DenseVector<Archs::CPU, DT_> c_local(size);
      DenseVector<Arch_, DT_> d(size);
      DenseVector<Archs::CPU, DT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i % 100 * DT_(1.234)));
        b_local(i, DT_(2 - i % 100));
        c_local(i, DT_(1 - i % 100));
        ref(i, c_local(i) * a_local(i) + b_local(i));
      }
      DenseVector<Arch_, DT_> a(a_local);
      DenseVector<Arch_, DT_> b(b_local);
      DenseVector<Arch_, DT_> c(c_local);

      Axpy<Arch_, BType_>::value(d, c, a, b);
      DenseVector<Archs::CPU, DT_> d_local(d);
      for (Index i(0) ; i < d_local.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(d_local(i), ref(i), 1e-2);
      Axpy<Arch_, BType_>::value(a, c, a, b);
      a_local = a;
      for (Index i(0) ; i < a_local.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(a_local(i), ref(i), 1e-2);
    }
  }
};
DVAxpyVTest<Archs::CPU, Archs::Generic, float> dv_axpy_v_test_float;
DVAxpyVTest<Archs::CPU, Archs::Generic, double> dv_axpy_v_test_double;
#ifdef FEAST_BACKENDS_CUDA
DVAxpyVTest<Archs::GPU, Archs::CUDA, float> gpu_dv_axpy_v_test_float;
DVAxpyVTest<Archs::GPU, Archs::CUDA, double> gpu_dv_axpy_v_test_double;
#endif
