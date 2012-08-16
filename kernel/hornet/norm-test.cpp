#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/hornet/dense_vector.hpp>
#include <kernel/hornet/norm.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class DVNorm2Test
  : public TaggedTest<Arch_, DT_>
{

public:

  DVNorm2Test()
    : TaggedTest<Arch_, DT_>("dv_norm2_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e5 ; size*=2)
    {
      DenseVector<Arch_, DT_> a(size);
      DT_ ref(0);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_((i * DT_(1.234) / size)));
        ref += a(i) * a(i);
      }
      ref = sqrt(ref);

      DT_ c = Norm2<Arch_, BType_>::value(a);
      TEST_CHECK_EQUAL(c, ref);
    }
  }
};
DVNorm2Test<Archs::CPU, Archs::Generic, float> dv_norm2_test_float;
DVNorm2Test<Archs::CPU, Archs::Generic, double> dv_norm2_test_double;
