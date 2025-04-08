// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/hash.hpp>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the Hash class.
 *
 * \test Tests the Hash class.
 *
 * \author Peter Zajac
 */
class HashTest :
  public TestSystem::UnitTest
{
public:
  HashTest() :
    TestSystem::UnitTest("HashTest")
  {
  }

  virtual ~HashTest()
  {
  }

  virtual void run() const override
  {
    TEST_CHECK_EQUAL(Hash::crc32(43u, "The quick brown fox jumps over the lazy dog"), 0x414FA339);
  }
} hash_test;
