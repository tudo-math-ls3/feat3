#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/scale.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class DVScaleTest
  : public TaggedTest<Arch_, DT_>
{

public:

  DVScaleTest()
    : TaggedTest<Arch_, DT_>("dv_scale_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e5 ; size*=2)
    {
      DT_ s(4.321);
      DenseVector<Arch_, DT_> a(size);
      DenseVector<Arch_, DT_> b(size);
      DenseVector<Arch_, DT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i * DT_(1.234)));
        ref(i, a(i) * s);
      }

      Scale<Arch_, BType_>::value(b, a, s);
      TEST_CHECK_EQUAL(b, ref);
      Scale<Arch_, BType_>::value(a, a, s);
      TEST_CHECK_EQUAL(a, ref);
    }
  }
};
DVScaleTest<Archs::CPU, Archs::Generic, float> dv_scale_test_float;
DVScaleTest<Archs::CPU, Archs::Generic, double> dv_scale_test_double;
