#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/hornet/dense_vector.hpp>
#include <kernel/hornet/difference.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class DVDifferenceTest
  : public TaggedTest<Arch_, DT_>
{

public:

  DVDifferenceTest()
    : TaggedTest<Arch_, DT_>("dv_difference_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e7 ; size*=2)
    {
      DenseVector<Arch_, DT_> a(size);
      DenseVector<Arch_, DT_> a2(size);
      DenseVector<Arch_, DT_> b(size);
      DenseVector<Arch_, DT_> c(size);
      DenseVector<Arch_, DT_> ref(size);
      DenseVector<Arch_, DT_> ref2(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i * DT_(1.234)));
        a2(i, DT_(i * DT_(1.234)));
        b(i, DT_(size*2 - i));
        ref(i, a(i) - b(i));
        ref2(i, b(i) - a(i));
      }

      Difference<Arch_, BType_>::value(c, a, b);
      TEST_CHECK_EQUAL(c, ref);
      Difference<Arch_, BType_>::value(c, b, a);
      TEST_CHECK_EQUAL(c, ref2);
      Difference<Arch_, BType_>::value(a, a, b);
      TEST_CHECK_EQUAL(a, ref);
      Difference<Arch_, BType_>::value(b, a2, b);
      TEST_CHECK_EQUAL(b, ref);
    }
  }
};
DVDifferenceTest<Archs::CPU, Archs::Generic, float> dv_difference_test_float;
DVDifferenceTest<Archs::CPU, Archs::Generic, double> dv_difference_test_double;
