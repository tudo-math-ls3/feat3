#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dot_product.hpp>
#include <kernel/lafem/norm.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVDotProductTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVDotProductTest()
    : TaggedTest<Arch_, DT_, Algo_>("dv_dot_product_test")
  {
  }

  virtual void run() const
  {
    const DT_ eps = Math::pow(Math::Limits<DT_>::epsilon(), DT_(0.8));

    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      const DT_ den(DT_(1) / DT_(size));
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i+1) * den);    // a[i] = (i+1) / n
        b_local(i, DT_(1) / DT_(i+1)); // b[i] = 1 / (i+1)
      }

      DenseVector<Arch_, DT_> a(a_local);
      DenseVector<Arch_, DT_> b(b_local);

      // a*b = 1
      DT_ ref(DT_(1));
      DT_ c = DotProduct<Algo_>::value(a, b);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = DotProduct<Algo_>::value(b, a);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = DotProduct<Algo_>::value(b, b);
      ref = Norm2<Algo_>::value(b);
      ref *= ref;
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
    }
  }
};
DVDotProductTest<Mem::Main, Algo::Generic, float> dv_dot_product_test_float;
DVDotProductTest<Mem::Main, Algo::Generic, double> dv_dot_product_test_double;
#ifdef FEAST_GMP
//DVDotProductTest<Mem::Main, Algo::Generic, mpf_class> dv_dot_product_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_MKL
DVDotProductTest<Mem::Main, Algo::MKL, float> mkl_dv_dot_product_test_float;
DVDotProductTest<Mem::Main, Algo::MKL, double> mkl_dv_dot_product_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVDotProductTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_dot_product_test_float;
DVDotProductTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_dot_product_test_double;
#endif
