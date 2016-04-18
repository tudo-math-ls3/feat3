#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/container_main_wrapper.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for the container main wrapper class.
 *
 * \test test description missing
 *
 * \tparam Mem_
 * description missing
 *
 * \tparam DT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename Mem_,
  typename DT_,
  typename IT_>
class ContainerMainWrapperTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  ContainerMainWrapperTest()
    : FullTaggedTest<Mem_, DT_, IT_>("ContainerMainWrapperTest")
  {
  }

  virtual void run() const override
  {
    DenseVector<Mem_, DT_, IT_> a(10, DT_(4711));
    {
      ContainerMainWrapper<decltype(a)> w(a);
      (*w).elements()[5] = DT_(42);
    }
    TEST_CHECK_EQUAL(a(4), DT_(4711));
    TEST_CHECK_EQUAL(a(5), DT_(42));
  }
};
ContainerMainWrapperTest<Mem::Main, float, unsigned int> cpu_container_main_wrapper_test_float_uint;
ContainerMainWrapperTest<Mem::Main, double, unsigned int> cpu_container_main_wrapper_test_double_uint;
ContainerMainWrapperTest<Mem::Main, float, unsigned long> cpu_container_main_wrapper_test_float_ulong;
ContainerMainWrapperTest<Mem::Main, double, unsigned long> cpu_container_main_wrapper_test_double_ulong;
#ifdef FEAST_HAVE_QUADMATH
ContainerMainWrapperTest<Mem::Main, __float128, unsigned int> cpu_container_main_wrapper_test_float128_uint;
ContainerMainWrapperTest<Mem::Main, __float128, unsigned long> cpu_container_main_wrapper_test_float128_ulong;
#endif
#ifdef FEAST_BACKENDS_CUDA
ContainerMainWrapperTest<Mem::CUDA, float, unsigned int> cuda_container_main_wrapper_test_float_uint;
ContainerMainWrapperTest<Mem::CUDA, double, unsigned int> cuda_container_main_wrapper_test_double_uint;
ContainerMainWrapperTest<Mem::CUDA, float, unsigned long> cuda_container_main_wrapper_test_float_ulong;
ContainerMainWrapperTest<Mem::CUDA, double, unsigned long> cuda_container_main_wrapper_test_double_ulong;
#endif
