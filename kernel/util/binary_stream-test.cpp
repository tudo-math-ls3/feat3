// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/util/binary_stream.hpp>
#include <stdint.h>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the BinaryStream class.
 *
 * \test Tests the BinaryStream class template.
 *
 * \author Peter Zajac
 */
class BinaryStreamTest :
  public TestSystem::UnitTest
{
public:
  BinaryStreamTest() :
    TestSystem::UnitTest("binary_stream_test")
  {
  }

  virtual ~BinaryStreamTest()
  {
  }

  static bool little_endian()
  {
    int32_t x = 1;
    void * xv = (void*)&x;
    int16_t * x16 = (int16_t*)xv;
    return *x16 == 1;
  }

  virtual void run() const override
  {
    BinaryStream bs;

    // write something to the stream
    const uint32_t db = 0xDEADBEEF;
    bs.write((const char*)&db, sizeof(uint32_t));

    // check stream position
    TEST_CHECK_EQUAL(size_t(bs.tellg()), sizeof(uint32_t));

    // seek back 2 bytes
    bs.seekp(-2, std::ios_base::cur);

    // read 2 bytes
    uint16_t d = 0;
    bs.read((char*)&d, sizeof(uint16_t));
    if(little_endian())
      TEST_CHECK_EQUAL(d, 0xDEAD); // high word
    else
      TEST_CHECK_EQUAL(d, 0xBEEF); // low word

    // write something else to the stream
    const uint32_t lc = 0x1337C0DE;
    bs.write((const char*)&lc, sizeof(uint32_t));

    // check stream position
    TEST_CHECK_EQUAL(size_t(bs.tellg()), 2*sizeof(uint32_t));

    // make a copy of the stream by writing
    BinaryStream bs_copy1;
    bs.write_stream(bs_copy1);

    // check for equal size
    TEST_CHECK_EQUAL(bs.size(), bs_copy1.size());

    // check for equal content
    for(std::size_t i(0); i < 2*sizeof(uint32_t); ++i)
      TEST_CHECK_EQUAL(bs.data()[i], bs_copy1.data()[i]);

    // make a copy of the stream by reading
    BinaryStream bs_copy2;
    bs_copy2.read_stream(bs);

    // check for equal size
    TEST_CHECK_EQUAL(bs.size(), bs_copy2.size());

    // check for equal content
    for(std::size_t i(0); i < 2*sizeof(uint32_t); ++i)
      TEST_CHECK_EQUAL(bs.data()[i], bs_copy2.data()[i]);

    // seek to begin
    bs.seekg(0);

    // read the whole stream
    uint32_t dblc[2] = {0, 0};
    bs.read((char*)dblc, 2*sizeof(uint32_t));
    TEST_CHECK_EQUAL(dblc[0], db);
    TEST_CHECK_EQUAL(dblc[1], lc);

    // check stream size
    TEST_CHECK_EQUAL(bs.size(), std::streamsize(8));

    // clear stream and check size
    bs.clear();
    TEST_CHECK_EQUAL(bs.size(), std::streamsize(0));
  }
} binary_stream_test;
