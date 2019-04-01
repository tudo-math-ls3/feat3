// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/util/memory_usage.hpp>
#include <thread>
#include <chrono>

using namespace FEAT;
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

  virtual ~MemoryUsageTest()
  {
  }

  virtual void run() const override
  {
    unsigned int * t = new unsigned int[123456];
    t[5] = 2;
    std::chrono::seconds sec(t[5]);
    std::this_thread::sleep_for(sec);

    MemoryUsage mu;
    TEST_CHECK_NOT_EQUAL(mu.get_current_physical(), 0);
    TEST_CHECK_NOT_EQUAL(mu.get_peak_physical(), 0);
    TEST_CHECK_NOT_EQUAL(mu.get_current_virtual(), 0);
    TEST_CHECK_NOT_EQUAL(mu.get_peak_virtual(), 0);
    TEST_CHECK_EQUAL(mu.get_current_swap(), 0);

    delete[] t;
  }
};
#if !defined(__linux) && defined(__unix)
// the memory usage on bsd systems is reported as zero in non interactive mode, but positive in interactive mode; thus it is useful but
// the test is disabled until somebody really uses a bsd system in his day to day work.
/// \todo fix memory usage report in non interactive mode on bsd systems
#else
MemoryUsageTest memory_usage_test;
#endif
