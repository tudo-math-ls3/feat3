// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/random.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Basic Memory::Arbiter test
 *
 * \author Peter Zajac
 */
class MemoryArbiterTest :
  public TestSystem::UnitTest
{
public:
  MemoryArbiterTest() :
    TestSystem::UnitTest("MemoryArbiterTest")
  {
  }

  void test_basic() const
  {
    static constexpr Index n = 17;

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));

    Memory::Arbiter a(n*sizeof(Index));

    TEST_CHECK_EQUAL(a.bytes(), std::size_t(n*sizeof(Index)));

    TEST_CHECK(a.location() == Memory::Location::none);
    TEST_CHECK_EQUAL(a.bytes(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(a.block_bytes(), std::size_t(n*sizeof(Index)));

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));

    {
      Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < n; ++i)
        v[i] = i;
    }

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));

    a.release();

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));
  }

  void test_format() const
  {
    static constexpr Index n = 17;

    Memory::Arbiter a(n*sizeof(Index), Memory::Location::main, Memory::Init::format_to_zero);

    TEST_CHECK(a.location() == Memory::Location::main);

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));

    {
      const Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < n; ++i)
      {
        TEST_CHECK_EQUAL(v(i), Index(0));
      }
    }

    a = Memory::Arbiter(n*sizeof(Index), Memory::Location::main, Memory::Init::format_to_one);

    TEST_CHECK(a.location() == Memory::Location::main);

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));

    {
      const Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < n; ++i)
      {
        TEST_CHECK_EQUAL(v(i), ~Index(0));
      }
    }

    a.format(); // to zero

    {
      const Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < n; ++i)
      {
        TEST_CHECK_EQUAL(v(i), Index(0));
      }
    }

    a.format(Index(0xDEADBEEF), Memory::Location::main);

    {
      const Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < n; ++i)
      {
        TEST_CHECK_EQUAL(v(i), Index(0xDEADBEEF));
      }
    }

    a.format_random();

    {
      const Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::read));
      Index all_or = Index(0);
      Index all_and = ~Index(0);
      for(Index i = 0; i < n; ++i)
      {
        all_or |= v(i);
        all_and &= v(i);
      }
      // if all OR-ed is 0, then all v(i) bits are zero and therefore not random
      TEST_CHECK_NOT_EQUAL(all_or, Index(0));
      // if all AND-ed is ~0, then all v(i) bits are one and therefore not random
      TEST_CHECK_NOT_EQUAL(all_and, ~Index(0));
    }
  }

  void test_invalid_attach() const
  {
    static constexpr std::size_t n = 27;

    Memory::Arbiter a(n, Memory::Location::main, Memory::Init::format_to_zero);

    // invalid offset
    TEST_CHECK_THROWS(a.attach(n + 1u), Memory::RangeException);

    // invalid size
    TEST_CHECK_THROWS(a.attach(n - 1u, std::size_t(3u)), Memory::RangeException);

    Memory::Arbiter b, c;

    TEST_CHECK_THROWS(b.attach(a, n + 1u), Memory::RangeException);

    TEST_CHECK_THROWS(c.attach(a, n - 1u, std::size_t(3u)), Memory::RangeException);
  }

  void test_range_attach() const
  {
    static constexpr Index n = 17;
    static constexpr Index off_b = 3;
    static constexpr Index len_b = 5;
    static constexpr Index off_c = 10;
    static constexpr Index len_c = 4;

    Memory::Arbiter a(n*sizeof(Index));

    TEST_CHECK_EQUAL(a.bytes(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(a.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < n; ++i)
        v[i] = n + i;
    }

    Memory::Arbiter b = a.attach(off_b*sizeof(Index), len_b*sizeof(Index));

    TEST_CHECK_EQUAL(b.bytes(), std::size_t(len_b*sizeof(Index)));
    TEST_CHECK_EQUAL(b.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      const Memory::TypedView<Index> v(b.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_b; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_b + i);
      }
    }

    {
      Memory::TypedView<Index> v(b.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < len_b; ++i)
      {
        v[i] = n - i;
      }
    }

    Memory::Arbiter c;
    c.attach(a, off_c * sizeof(Index), len_c*sizeof(Index));

    TEST_CHECK_EQUAL(c.bytes(), std::size_t(len_c*sizeof(Index)));
    TEST_CHECK_EQUAL(c.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      const Memory::TypedView<Index> v(c.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c + i);
      }
    }

    {
      Memory::TypedView<Index> v(c.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < len_c; ++i)
      {
        v[i] = 3u*n - i;
      }
    }

    {
      const Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::read));
      // [0, off_b): original values of a
      for(Index i = 0; i < off_b; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + i);
      }
      // [off_b, off_b+len_b): overwritten by b
      for(Index i = 0; i < len_b; ++i)
      {
        TEST_CHECK_EQUAL(v(off_b + i), n - i);
      }
      // [off_b+len_b,off_c): original values of a
      for(Index i = off_b + len_b; i < off_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + i);
      }
      // [off_c, off_c+len_c): overwritten by c
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(off_c + i), 3u*n - i);
      }
      // [off_c+len_c,n): original values of a
      for(Index i = off_c + len_c; i < n; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + i);
      }
    }
  }

  void test_subrange_attach() const
  {
    static constexpr Index n = 17;
    static constexpr Index off_b = 3;
    static constexpr Index len_b = 12;
    static constexpr Index off_c_b = 2;
    static constexpr Index len_c = 7;
    static constexpr Index off_c_a = off_b + off_c_b;
    static_assert(off_c_a + len_c <= n);
    static_assert(off_c_b + len_c <= len_b);

    Memory::Arbiter a(n*sizeof(Index));

    TEST_CHECK_EQUAL(a.bytes(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(a.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < n; ++i)
        v[i] = n + i;
    }

    Memory::Arbiter b = a.attach(off_b*sizeof(Index), len_b*sizeof(Index));

    TEST_CHECK_EQUAL(b.bytes(), std::size_t(len_b*sizeof(Index)));
    TEST_CHECK_EQUAL(b.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      const Memory::TypedView<Index> v(b.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_b; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_b + i);
      }
    }

    Memory::Arbiter c1 = b.attach(off_c_b*sizeof(Index), len_c*sizeof(Index));

    TEST_CHECK_EQUAL(c1.bytes(), std::size_t(len_c*sizeof(Index)));
    TEST_CHECK_EQUAL(c1.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      const Memory::TypedView<Index> v(c1.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c_a + i);
      }
    }

    Memory::Arbiter c2;
    c2.attach(b, off_c_b*sizeof(Index), len_c*sizeof(Index));

    TEST_CHECK_EQUAL(c2.bytes(), std::size_t(len_c*sizeof(Index)));
    TEST_CHECK_EQUAL(c2.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      const Memory::TypedView<Index> v(c2.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c_a + i);
      }
    }

    Memory::Arbiter c3  = b.attach();
    c3.attach(c3, off_c_b*sizeof(Index), len_c*sizeof(Index)); // self-attach

    TEST_CHECK_EQUAL(c3.bytes(), std::size_t(len_c*sizeof(Index)));
    TEST_CHECK_EQUAL(c3.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      const Memory::TypedView<Index> v(c3.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c_a + i);
      }
    }
  }

  void test_invalid_clone() const
  {
    static constexpr std::size_t n = 27;

    Memory::Arbiter a(n, Memory::Location::main, Memory::Init::format_to_zero);

    // invalid offset
    TEST_CHECK_THROWS(a.clone(n + 1u), Memory::RangeException);

    // invalid size
    TEST_CHECK_THROWS(a.clone(n - 1u, std::size_t(3u)), Memory::RangeException);

    Memory::Arbiter b, c;

    TEST_CHECK_THROWS(b.clone(a, n + 1u), Memory::RangeException);

    TEST_CHECK_THROWS(c.clone(a, n - 1u, std::size_t(3u)), Memory::RangeException);
  }

  void test_range_clone() const
  {
    static constexpr Index n = 17;
    static constexpr Index off_b = 3;
    static constexpr Index len_b = 5;
    static constexpr Index off_c = 10;
    static constexpr Index len_c = 4;

    Memory::Arbiter a(n*sizeof(Index));

    TEST_CHECK_EQUAL(a.bytes(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(a.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < n; ++i)
        v[i] = n + i;
    }

    // clone entire block
    Memory::Arbiter b1 = a.clone(off_b*sizeof(Index), len_b*sizeof(Index), true);

    TEST_CHECK_EQUAL(b1.bytes(), std::size_t(len_b*sizeof(Index)));
    TEST_CHECK_EQUAL(b1.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      const Memory::TypedView<Index> v(b1.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_b; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_b + i);
      }
    }

    {
      Memory::TypedView<Index> v(b1.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < len_b; ++i)
      {
        v[i] = n - i;
      }
    }

    // clone range block
    Memory::Arbiter b2 = a.clone(off_b*sizeof(Index), len_b*sizeof(Index), false);

    TEST_CHECK_EQUAL(b2.bytes(), std::size_t(len_b*sizeof(Index)));
    TEST_CHECK_EQUAL(b2.block_bytes(), std::size_t(len_b*sizeof(Index)));

    {
      const Memory::TypedView<Index> v(b2.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_b; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_b + i);
      }
    }

    {
      Memory::TypedView<Index> v(b2.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < len_b; ++i)
      {
        v[i] = n - i;
      }
    }

    // clone entire block
    Memory::Arbiter c1;
    c1.clone(a, off_c * sizeof(Index), len_c*sizeof(Index), true);

    TEST_CHECK_EQUAL(c1.bytes(), std::size_t(len_c*sizeof(Index)));
    TEST_CHECK_EQUAL(c1.block_bytes(), std::size_t(n*sizeof(Index)));

    {
      const Memory::TypedView<Index> v(c1.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c + i);
      }
    }

    {
      Memory::TypedView<Index> v(c1.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < len_c; ++i)
      {
        v[i] = 3u*n - i;
      }
    }

    // clone range block
    Memory::Arbiter c2;
    c2.clone(a, off_c * sizeof(Index), len_c*sizeof(Index), false);

    TEST_CHECK_EQUAL(c2.bytes(), std::size_t(len_c*sizeof(Index)));
    TEST_CHECK_EQUAL(c2.block_bytes(), std::size_t(len_c*sizeof(Index)));

    {
      const Memory::TypedView<Index> v(c2.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c + i);
      }
    }

    {
      Memory::TypedView<Index> v(c2.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < len_c; ++i)
      {
        v[i] = 3u*n - i;
      }
    }

    // a should be unmodified
    {
      const Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < n; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + i);
      }
    }
  }

  void test_subrange_clone_full() const
  {
    static constexpr Index n = 17;
    static constexpr Index off_b = 3;
    static constexpr Index len_b = 12;
    static constexpr Index off_c_b = 2;
    static constexpr Index len_c = 7;
    static constexpr Index off_c_a = off_b + off_c_b;
    static_assert(off_c_a + len_c <= n);
    static_assert(off_c_b + len_c <= len_b);

    Memory::Arbiter a(n*sizeof(Index));

    TEST_CHECK_EQUAL(a.bytes(), n*sizeof(Index));
    TEST_CHECK_EQUAL(a.block_bytes(), n*sizeof(Index));

    {
      Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < n; ++i)
        v[i] = n + i;
    }

    Memory::Arbiter b = a.clone(off_b*sizeof(Index), len_b*sizeof(Index), true);

    TEST_CHECK_EQUAL(b.bytes(), len_b*sizeof(Index));
    TEST_CHECK_EQUAL(b.block_bytes(), n*sizeof(Index));

    {
      const Memory::TypedView<Index> v(b.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_b; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_b + i);
      }
    }

    Memory::Arbiter c1 = b.clone(off_c_b*sizeof(Index), len_c*sizeof(Index), true);

    TEST_CHECK_EQUAL(c1.bytes(), len_c*sizeof(Index));
    TEST_CHECK_EQUAL(c1.block_bytes(), n*sizeof(Index));

    {
      const Memory::TypedView<Index> v(c1.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c_a + i);
      }
    }

    Memory::Arbiter c2;
    c2.clone(b, off_c_b*sizeof(Index), len_c*sizeof(Index), true);

    TEST_CHECK_EQUAL(c2.bytes(), len_c*sizeof(Index));
    TEST_CHECK_EQUAL(c2.block_bytes(), n*sizeof(Index));

    {
      const Memory::TypedView<Index> v(c2.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c_a + i);
      }
    }

    Memory::Arbiter c3  = b.clone(std::size_t(0), Memory::Arbiter::full_size, true);

    {
      const Memory::TypedView<Index> v(c3.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_b; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_b + i);
      }
    }

    c3.clone(c3, off_c_b*sizeof(Index), len_c*sizeof(Index), true); // self-clone

    TEST_CHECK_EQUAL(c3.bytes(), len_c*sizeof(Index));
    TEST_CHECK_EQUAL(c3.block_bytes(), n*sizeof(Index));

    {
      const Memory::TypedView<Index> v(c3.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c_a + i);
      }
    }
  }

  void test_subrange_clone_range() const
  {
    static constexpr Index n = 17;
    static constexpr Index off_b = 3;
    static constexpr Index len_b = 12;
    static constexpr Index off_c_b = 2;
    static constexpr Index len_c = 7;
    static constexpr Index off_c_a = off_b + off_c_b;
    static_assert(off_c_a + len_c <= n);
    static_assert(off_c_b + len_c <= len_b);

    Memory::Arbiter a(n*sizeof(Index));

    TEST_CHECK_EQUAL(a.bytes(), n*sizeof(Index));
    TEST_CHECK_EQUAL(a.block_bytes(), n*sizeof(Index));

    {
      Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < n; ++i)
        v[i] = n + i;
    }

    Memory::Arbiter b = a.clone(off_b*sizeof(Index), len_b*sizeof(Index), false);

    TEST_CHECK_EQUAL(b.bytes(), len_b*sizeof(Index));
    TEST_CHECK_EQUAL(b.block_bytes(), len_b*sizeof(Index));

    {
      const Memory::TypedView<Index> v(b.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_b; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_b + i);
      }
    }

    Memory::Arbiter c1 = b.clone(off_c_b*sizeof(Index), len_c*sizeof(Index), false);

    TEST_CHECK_EQUAL(c1.bytes(), len_c*sizeof(Index));
    TEST_CHECK_EQUAL(c1.block_bytes(), len_c*sizeof(Index));

    {
      const Memory::TypedView<Index> v(c1.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c_a + i);
      }
    }

    Memory::Arbiter c2;
    c2.clone(b, off_c_b*sizeof(Index), len_c*sizeof(Index), false);

    TEST_CHECK_EQUAL(c2.bytes(), len_c*sizeof(Index));
    TEST_CHECK_EQUAL(c2.block_bytes(), len_c*sizeof(Index));

    {
      const Memory::TypedView<Index> v(c2.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c_a + i);
      }
    }

    Memory::Arbiter c3  = b.clone(std::size_t(0), Memory::Arbiter::full_size, false);

    {
      const Memory::TypedView<Index> v(c3.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_b; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_b + i);
      }
    }

    c3.clone(c3, off_c_b*sizeof(Index), len_c*sizeof(Index), false); // self-clone

    TEST_CHECK_EQUAL(c3.bytes(), len_c*sizeof(Index));
    TEST_CHECK_EQUAL(c3.block_bytes(), len_c*sizeof(Index));

    {
      const Memory::TypedView<Index> v(c3.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < len_c; ++i)
      {
        TEST_CHECK_EQUAL(v(i), n + off_c_a + i);
      }
    }
  }

  void test_read_write_overlap() const
  {
    // initialize in cuda with manual memcopy, create write view in main and derive read view in main
    static constexpr Index n = 17;

    Memory::Arbiter a(n*sizeof(Index), Memory::Location::main, Memory::Init::format_to_zero);

    {
      Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<Index> w(a.view(Memory::Location::main, Memory::Access::read));
      TEST_CHECK_EQUAL(w.raw_r(), v.raw_r());
    }

    a.format();

    Memory::Arbiter b(n*sizeof(Index), Memory::Location::main, Memory::Init::format_to_one);

    {
      Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<Index> w(b.view(Memory::Location::main, Memory::Access::read));
      TEST_CHECK_NOT_EQUAL(w.raw_r(), v.raw_r());
    }
  }

  void test_write_overlap() const
  {
    static constexpr Index n = 17;

    Memory::Arbiter a(n*sizeof(Index));

    {
      Memory::TypedView<Index> v1(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<Index> v2(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      for(Index i = 0; i < n; i += 2)
        v1[i] = i;
      for(Index i = 1; i < n; i += 2)
        v2[i] = i;
    }

    {
      const Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < n; ++i)
      {
        TEST_CHECK_EQUAL(v(i), i);
      }
    }
  }

  void test_access_exceptions() const
  {
    static constexpr Index n = 17;

    Memory::Arbiter a(n*sizeof(Index), Memory::Location::main, Memory::Init::format_to_zero);

    // read in read (ok)
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::read));
      Memory::View w(a.view(Memory::Location::main, Memory::Access::read));
      TEST_CHECK_EQUAL(v.raw_r(), w.raw_r());
    }

    // read in write
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write));
      TEST_CHECK_THROWS(a.view(Memory::Location::main, Memory::Access::read), Memory::AccessException);
    }

    // write in read
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::read));
      TEST_CHECK_THROWS(a.view(Memory::Location::main, Memory::Access::write), Memory::AccessException);
    }

    // write in write
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write));
      TEST_CHECK_THROWS(a.view(Memory::Location::main, Memory::Access::write), Memory::AccessException);
    }

    // overlap write in write
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write));
      TEST_CHECK_THROWS(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap), Memory::AccessException);
    }

    // write in overlap write
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      TEST_CHECK_THROWS(a.view(Memory::Location::main, Memory::Access::write), Memory::AccessException);
    }

    // overlap write in overlap write (ok)
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::View w(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      TEST_CHECK_EQUAL(v.raw_w(), w.raw_w());
    }

    // read in overlap write (ok)
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::View w(a.view(Memory::Location::main, Memory::Access::read));
      TEST_CHECK_EQUAL(v.raw_r(), w.raw_r());
    }

    // overlap write in read
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::read));
      TEST_CHECK_THROWS(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap), Memory::AccessException);
    }
  }

  virtual void run() const override
  {
    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));

    test_basic();
    test_format();
    test_invalid_attach();
    test_range_attach();
    test_subrange_attach();
    test_invalid_clone();
    test_range_clone();
    test_subrange_clone_full();
    test_subrange_clone_range();
    test_read_write_overlap();
    test_write_overlap();
    test_access_exceptions();

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));
  }
} memory_arbiter_test;

#ifdef FEAT_HAVE_CUDA

/**
 * \brief Memory::Arbiter CUDA test
 *
 * \author Peter Zajac
 */
class MemoryArbiterCUDATest :
  public TestSystem::UnitTest
{
public:
  MemoryArbiterCUDATest() :
    TestSystem::UnitTest("MemoryArbiterCUDATest")
  {
  }

  void test_basic() const
  {
    static constexpr Index n = 17;

    Memory::Arbiter a(n*sizeof(Index));

    TEST_CHECK(a.location() == Memory::Location::none);

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));

    {
      Memory::TypedView<Index> view(a.view(Memory::Location::main, Memory::Access::write));
    }

    TEST_CHECK(a.location() == Memory::Location::main);

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));

    {
      Memory::TypedView<Index> view(a.view(Memory::Location::cuda, Memory::Access::read));
    }

    TEST_CHECK(a.location() == (Memory::Location::main|Memory::Location::cuda));

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(n*sizeof(Index)));

    Memory::Manager::free_all_unmapped(Memory::Location::cuda);

    TEST_CHECK(a.location() == Memory::Location::main);

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));

    {
      Memory::TypedView<Index> view(a.view(Memory::Location::cuda, Memory::Access::write));
    }

    TEST_CHECK(a.location() == Memory::Location::cuda);

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(n*sizeof(Index)));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(n*sizeof(Index)));
  }

  void test_main_to_cuda() const
  {
    // initialize in main, create read view in cuda, manual memcopy to main
    static constexpr Index n = 17;
    Index x[n];

    Memory::Arbiter a(n*sizeof(Index));

    {
      Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < n; ++i)
        v[i] = i;
    }
    {
      Memory::TypedView<Index> r(a.view(Memory::Location::cuda, Memory::Access::read));
      Memory::memcopy_cuda_to_main(x, r.get_r(), n*sizeof(Index));
    }
    for(Index i = 0; i < n; ++i)
    {
      TEST_CHECK_EQUAL(x[i], i);
    }
  }

  void test_cuda_to_main() const
  {
    // initialize in cuda with manual memcopy, create read view in main
    static constexpr Index n = 17;
    Index x[n];
    Memory::Arbiter a(n*sizeof(Index));
    {
      Memory::TypedView<Index> v(a.view(Memory::Location::cuda, Memory::Access::write));
      for(Index i = 0; i < n; ++i)
        x[i] = i;
      Memory::memcopy_main_to_cuda(v.raw_w(), x, n*sizeof(Index));
    }
    {
      Memory::TypedView<Index> r(a.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < n; ++i)
      {
        TEST_CHECK_EQUAL(r(i), i);
      }
    }
  }

  void test_derive_view_main() const
  {
    // initialize in cuda with manual memcopy, create write view in main and derive read view in main
    static constexpr Index n = 17;
    Index x[n];

    Memory::Arbiter a(n*sizeof(Index));
    {
      Memory::TypedView<Index> v(a.view(Memory::Location::cuda, Memory::Access::write));
      for(Index i = 0; i < n; ++i)
        x[i] = i;
      Memory::memcopy_main_to_cuda(v.raw_w(), x, n*sizeof(Index));
    }
    {
      Memory::TypedView<Index> w(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<Index> r(a.view(Memory::Location::main, Memory::Access::read));
      for(Index i = 0; i < n; ++i)
      {
        TEST_CHECK_EQUAL(r(i), i);
      }
    }
  }

  void test_derive_view_cuda() const
  {
    // initialize in main, create write view in cuda, derive read view in cuda, manual memcopy to main
    static constexpr Index n = 17;
    Index x[n];

    Memory::Arbiter a(n*sizeof(Index));
    {
      Memory::TypedView<Index> v(a.view(Memory::Location::main, Memory::Access::write));
      for(Index i = 0; i < n; ++i)
        v[i] = i;
    }
    {
      Memory::TypedView<Index> w(a.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<Index> r(a.view(Memory::Location::cuda, Memory::Access::read));
      Memory::memcopy_cuda_to_main(x, r.get_r(), n*sizeof(Index));
    }
    for(Index i = 0; i < n; ++i)
    {
      TEST_CHECK_EQUAL(x[i], i);
    }
  }

  void test_access_exceptions() const
  {
    static constexpr Index n = 17;

    Memory::Arbiter a(n*sizeof(Index), Memory::Location::main, Memory::Init::format_to_zero);

    // read in read (ok)
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::read));
      Memory::View w(a.view(Memory::Location::cuda, Memory::Access::read));
      TEST_CHECK_NOT_EQUAL(v.raw_r(), w.raw_r());
    }

    // read in write
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write));
      TEST_CHECK_THROWS(a.view(Memory::Location::cuda, Memory::Access::read), Memory::AccessException);
    }

    // write in read
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::read));
      TEST_CHECK_THROWS(a.view(Memory::Location::cuda, Memory::Access::write), Memory::AccessException);
    }

    // write in write
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write));
      TEST_CHECK_THROWS(a.view(Memory::Location::cuda, Memory::Access::write), Memory::AccessException);
    }

    // overlap write in write
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write));
      TEST_CHECK_THROWS(a.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap), Memory::AccessException);
    }

    // write in overlap write
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      TEST_CHECK_THROWS(a.view(Memory::Location::cuda, Memory::Access::write), Memory::AccessException);
    }

    // overlap write in overlap write (different locations)
    {
      Memory::View v(a.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      TEST_CHECK_THROWS(a.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap), Memory::AccessException);
    }

    // read in overlap write (ok)
    {
      Memory::View v(a.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      Memory::View w(a.view(Memory::Location::cuda, Memory::Access::read));
      TEST_CHECK_EQUAL(v.raw_r(), w.raw_r());
    }

    // overlap write in read
    {
      Memory::View v(a.view(Memory::Location::cuda, Memory::Access::read));
      TEST_CHECK_THROWS(a.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap), Memory::AccessException);
    }
  }

  virtual void run() const override
  {
    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));

    test_basic();
    test_main_to_cuda();
    test_cuda_to_main();
    test_derive_view_main();
    test_derive_view_cuda();
    test_access_exceptions();

    TEST_CHECK_EQUAL(Memory::Manager::bytes_reserved(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_main(), std::size_t(0));
    TEST_CHECK_EQUAL(Memory::Manager::bytes_allocated_cuda(), std::size_t(0));
  }
} memory_arbiter_cuda_test;

#endif // FEAT_HAVE_CUDA
