#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/element_product.hpp>
#include <kernel/lafem/algorithm.hpp>

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
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> ref2(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i * DT_(1.234)));
        b_local(i, DT_(size*2 - i));
        ref(i, a_local(i) * b_local(i));
        ref2(i, a_local(i) * a_local(i));
      }

      DenseVector<Arch_, DT_> a(size);
      copy(a, a_local);
      DenseVector<Arch_, DT_> b(size);
      copy(b, b_local);
      DenseVector<Arch_, DT_> c(size);

      ElementProduct<Arch_, BType_>::value(c, a, b);
      copy(result_local, c);
      TEST_CHECK_EQUAL(result_local, ref);

      ElementProduct<Arch_, BType_>::value(a, a, b);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(a, a_local);
      ElementProduct<Arch_, BType_>::value(b, a, b);
      copy(result_local, b);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(b, b_local);
      ElementProduct<Arch_, BType_>::value(a, a, a);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref2);
    }
  }
};
DVElementProductTest<Mem::Main, Algo::Generic, float> dv_element_product_test_float;
DVElementProductTest<Mem::Main, Algo::Generic, double> dv_element_product_test_double;
#ifdef FEAST_BACKENDS_CUDA
DVElementProductTest<Mem::CUDA, Algo::CUDA, float> gpu_dv_element_product_test_float;
DVElementProductTest<Mem::CUDA, Algo::CUDA, double> gpu_dv_element_product_test_double;
#endif
