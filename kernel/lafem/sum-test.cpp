#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sum.hpp>
#include <kernel/lafem/algorithm.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class DVSumTest
  : public TaggedTest<Arch_, DT_>
{

public:

  DVSumTest()
    : TaggedTest<Arch_, DT_>("dv_sum_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e5 ; size*=2)
    {
      DenseVector<Archs::CPU, DT_> a_local(size);
      DenseVector<Archs::CPU, DT_> b_local(size);
      DenseVector<Archs::CPU, DT_> ref(size);
      DenseVector<Archs::CPU, DT_> ref2(size);
      DenseVector<Archs::CPU, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i * DT_(1.234)));
        b_local(i, DT_(size*2 - i));
        ref(i, a_local(i) + b_local(i));
        ref2(i, a_local(i) + a_local(i));
      }

      DenseVector<Arch_, DT_> a(size);
      copy(a, a_local);
      DenseVector<Arch_, DT_> b(size);
      copy(b, b_local);
      DenseVector<Arch_, DT_> c(size);

      Sum<Arch_, BType_>::value(c, a, b);
      copy(result_local, c);
      TEST_CHECK_EQUAL(result_local, ref);

      Sum<Arch_, BType_>::value(a, a, b);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(a, a_local);
      Sum<Arch_, BType_>::value(b, a, b);
      copy(result_local, b);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(b, b_local);
      Sum<Arch_, BType_>::value(a, a, a);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref2);
    }
  }
};
DVSumTest<Archs::CPU, Archs::Generic, float> dv_sum_test_float;
DVSumTest<Archs::CPU, Archs::Generic, double> dv_sum_test_double;
