#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dot_product.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class DVDotProductTest
  : public TaggedTest<Arch_, DT_>
{

public:

  DVDotProductTest()
    : TaggedTest<Arch_, DT_>("dv_dot_product_test")
  {
  }

  virtual void run() const
  {
    const DT_ eps = std::pow(std::numeric_limits<DT_>::epsilon(), DT_(0.8));

    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Arch_, DT_> a(size);
      DenseVector<Arch_, DT_> b(size);
      const DT_ den(DT_(1) / DT_(size));
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i+1) * den);    // a[i] = (i+1) / n
        b(i, DT_(1) / DT_(i+1)); // b[i] = 1 / (i+1)
      }

      // a*b = 1
      const DT_ ref(DT_(1));
      DT_ c = DotProduct<Arch_, BType_>::value(a, b);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = DotProduct<Arch_, BType_>::value(b, a);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
    }
  }
};
DVDotProductTest<Archs::CPU, Archs::Generic, float> dv_dot_product_test_float;
DVDotProductTest<Archs::CPU, Archs::Generic, double> dv_dot_product_test_double;
