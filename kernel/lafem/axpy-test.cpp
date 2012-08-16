#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/axpy.hpp>

using namespace FEAST;
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
    DT_ s(4711.1);
    for (Index size(1) ; size < 1e5 ; size*=2)
    {
      DenseVector<Arch_, DT_> a(size);
      DenseVector<Arch_, DT_> b(size);
      DenseVector<Arch_, DT_> c(size);
      DenseVector<Arch_, DT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i * DT_(1.234)));
        b(i, DT_(size*2 - i));
        ref(i, s * a(i) + b(i));
      }

      Axpy<Arch_, BType_>::value(c, s, a, b);
      TEST_CHECK_EQUAL(c, ref);
      Axpy<Arch_, BType_>::value(a, s, a, b);
      TEST_CHECK_EQUAL(a, ref);
    }
  }
};
DVAxpyTest<Archs::CPU, Archs::Generic, float> dv_axpy_test_float;
DVAxpyTest<Archs::CPU, Archs::Generic, double> dv_axpy_test_double;

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
      DenseVector<Arch_, DT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i * DT_(1.234)));
        b(i, DT_(size*2 - i));
        c(i, DT_(size - i));
        ref(i, c(i) * a(i) + b(i));
      }

      Axpy<Arch_, BType_>::value(d, c, a, b);
      TEST_CHECK_EQUAL(d, ref);
      Axpy<Arch_, BType_>::value(a, c, a, b);
      TEST_CHECK_EQUAL(a, ref);
    }
  }
};
DVAxpyVTest<Archs::CPU, Archs::Generic, float> dv_axpy_v_test_float;
DVAxpyVTest<Archs::CPU, Archs::Generic, double> dv_axpy_v_test_double;
