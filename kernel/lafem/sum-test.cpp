#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sum.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class DVSumTest
  : public TaggedTest<Arch_, DT_>
{

public:

  DVSumTest()
    : TaggedTest<Arch_, DT_>("dv_sum_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e5 ; size*=2)
    {
      DenseVector<Arch_, DT_> a(size);
      DenseVector<Arch_, DT_> a2(size);
      DenseVector<Arch_, DT_> b(size);
      DenseVector<Arch_, DT_> c(size);
      DenseVector<Arch_, DT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i * DT_(1.234)));
        a2(i, DT_(i * DT_(1.234)));
        b(i, DT_(size*2 - i));
        ref(i, a(i) + b(i));
      }

      Sum<Arch_, BType_>::value(c, b, a);
      TEST_CHECK_EQUAL(c, ref);
      Sum<Arch_, BType_>::value(a, a, b);
      TEST_CHECK_EQUAL(a, ref);
      Sum<Arch_, BType_>::value(b, a2, b);
      TEST_CHECK_EQUAL(b, ref);
    }
  }
};
DVSumTest<Archs::CPU, Archs::Generic, float> dv_sum_test_float;
DVSumTest<Archs::CPU, Archs::Generic, double> dv_sum_test_double;
