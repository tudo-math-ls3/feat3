#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/element_product.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class DVElementProductTest
  : public TaggedTest<Arch_, DT_>
{

public:

  DVElementProductTest()
    : TaggedTest<Arch_, DT_>("dv_element_product_test")
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
        ref(i, a(i) * b(i));
      }

      ElementProduct<Arch_, BType_>::value(c, b, a);
      TEST_CHECK_EQUAL(c, ref);
      ElementProduct<Arch_, BType_>::value(a, a, b);
      TEST_CHECK_EQUAL(a, ref);
      ElementProduct<Arch_, BType_>::value(b, a2, b);
      TEST_CHECK_EQUAL(b, ref);
    }
  }
};
DVElementProductTest<Archs::CPU, Archs::Generic, float> dv_element_product_test_float;
DVElementProductTest<Archs::CPU, Archs::Generic, double> dv_element_product_test_double;
