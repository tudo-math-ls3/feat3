#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/scale.hpp>
#include <kernel/lafem/algorithm.hpp>

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
      DT_ s(DT_(4.321));
      DenseVector<Archs::CPU, DT_> a_local(size);
      DenseVector<Archs::CPU, DT_> ref(size);
      DenseVector<Archs::CPU, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i * DT_(1.234)));
        ref(i, a_local(i) * s);
      }

      DenseVector<Arch_, DT_> a(size);
      copy(a, a_local);
      DenseVector<Arch_, DT_> b(size);

      Scale<Arch_, BType_>::value(b, a, s);
      copy(result_local, b);
      TEST_CHECK_EQUAL(result_local, ref);

      Scale<Arch_, BType_>::value(a, a, s);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref);
    }
  }
};
DVScaleTest<Archs::CPU, Archs::Generic, float> dv_scale_test_float;
DVScaleTest<Archs::CPU, Archs::Generic, double> dv_scale_test_double;
#ifdef FEAST_BACKENDS_CUDA
DVScaleTest<Archs::GPU, Archs::CUDA, float> cuda_dv_scale_test_float;
DVScaleTest<Archs::GPU, Archs::CUDA, double> cuda_dv_scale_test_double;
#endif
