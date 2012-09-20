#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/norm.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
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
    const DT_ eps = std::pow(std::numeric_limits<DT_>::epsilon(), DT_(0.8));

    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Archs::CPU, DT_> a_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        // a[i] = 1/sqrt(2^i) = (1/2)^(i/2)
        a_local(i, std::pow(DT_(0.5), DT_(0.5) * DT_(i)));
      }

      // ||a||_2 = sqrt(2 - 2^{1-n})
      const DT_ ref(std::sqrt(DT_(2) - std::pow(DT_(0.5), DT_(size-1))));

      DenseVector<Arch_, DT_> a(a_local);
      DT_ c = Norm2<Arch_, BType_>::value(a);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
    }
  }
};
DVNorm2Test<Archs::CPU, Archs::Generic, float> dv_norm2_test_float;
DVNorm2Test<Archs::CPU, Archs::Generic, double> dv_norm2_test_double;
#ifdef FEAST_BACKENDS_CUDA
DVNorm2Test<Archs::GPU, Archs::CUDA, float> cuda_dv_norm2_test_float;
DVNorm2Test<Archs::GPU, Archs::CUDA, double> cuda_dv_norm2_test_double;
#endif
