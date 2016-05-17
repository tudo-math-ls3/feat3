#include <test_system/test_system.hpp>
#include <kernel/util/memory_usage.hpp>

using namespace FEAST;
using namespace FEAST::Util;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for the memory usage method.
 *
 * \author Dirk Ribbrock
 */
class MemoryUsageTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  MemoryUsageTest() :
    TaggedTest<Archs::None, Archs::None>("memory_usage_test")
  {
  }

  virtual void run() const override
  {
    auto m = get_memory_usage();
    TEST_CHECK_NOT_EQUAL(m.current_physical, 0);
    TEST_CHECK_NOT_EQUAL(m.peak_physical, 0);
    TEST_CHECK_NOT_EQUAL(m.current_virtual, 0);
    TEST_CHECK_NOT_EQUAL(m.peak_virtual, 0);
    TEST_CHECK_EQUAL(m.current_swap, 0);
  }
} memory_usage_test;
