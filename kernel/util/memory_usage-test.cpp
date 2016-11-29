#include <test_system/test_system.hpp>
#include <kernel/util/memory_usage.hpp>
#include <unistd.h>

using namespace FEAT;
using namespace FEAT::Util;
using namespace FEAT::TestSystem;

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
    unsigned int * t = new unsigned int[123456];
    t[5] = 2;
    sleep(t[5]);

    auto m = get_memory_usage();
    TEST_CHECK_NOT_EQUAL(m.current_physical, 0);
    TEST_CHECK_NOT_EQUAL(m.peak_physical, 0);
    TEST_CHECK_NOT_EQUAL(m.current_virtual, 0);
    TEST_CHECK_NOT_EQUAL(m.peak_virtual, 0);
    TEST_CHECK_EQUAL(m.current_swap, 0);

    delete[] t;
  }
} memory_usage_test;
