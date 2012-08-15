#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/hornet/dense_vector.hpp>
#include <kernel/hornet/dot_product.hpp>

using namespace FEAST;
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
    for (Index size(1) ; size < 1e7 ; size*=2)
    {
      DenseVector<Arch_, DT_> a(size);
      DenseVector<Arch_, DT_> b(size);
      DT_ ref(0);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_((i * DT_(1.234)) / size));
        b(i, DT_((size*2 - i) / size));
        ref += a(i) * b(i);
      }

      DT_ c = DotProduct<Arch_, BType_>::value(a, b);
      TEST_CHECK_EQUAL(c, ref);
      c = DotProduct<Arch_, BType_>::value(b, a);
      TEST_CHECK_EQUAL(c, ref);
    }
  }
};
DVDotProductTest<Archs::CPU, Archs::Generic, float> dv_dot_product_test_float;
DVDotProductTest<Archs::CPU, Archs::Generic, double> dv_dot_product_test_double;
