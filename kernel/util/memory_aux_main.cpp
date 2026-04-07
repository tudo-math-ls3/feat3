// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/random.hpp>

// includes, system
#include <memory>
#include <cstring>

namespace FEAT
{
  namespace Memory
  {
    void* alloc_main(std::size_t bytes)
    {
      return malloc(bytes);
    }

    void free_main(void* ptr)
    {
      free(ptr);
    }

    void memset_main(void* ptr, int value, std::size_t bytes)
    {
      if(bytes <= std::size_t(0))
        return;

      XASSERT(ptr != nullptr);

      ::memset(ptr, value, bytes);
    }

    void memset_random_main(void* ptr, std::size_t bytes, std::uint64_t seed)
    {
      if(bytes <= std::size_t(0))
        return;

      XASSERT(ptr != nullptr);

      // zero seed? replace by pointer address then
      if(seed == 0u)
        seed = reinterpret_cast<std::uint64_t>(ptr);

      // initialize RNG
      Random rng(seed);

      // cast array to 64 bit int and fill them with random values
      std::uint64_t* x = reinterpret_cast<std::uint64_t*>(ptr);
      std::size_t n = bytes / 8u;
      for(std::size_t i = 0; i < n; ++i, ++x)
        *x = rng.next();

      // there may be a few (<8) bytes left, so format them via a memcpy
      std::size_t bytes_left = bytes & 0x7u;
      if(bytes_left > 0u)
      {
        std::uint64_t t = rng.next();
        ::memcpy(x, &t, bytes_left);
      }
    }

    void memcopy_main(void* dst, const void* src, std::size_t bytes)
    {
      if(bytes <= std::size_t(0))
        return;

      XASSERT(dst != nullptr);
      XASSERT(src != nullptr);

      if(dst == src)
        return;

      // check for overlapping memory regions
      XASSERTM((dst < src) || ((const char*)src + bytes <= dst), "source and destination memory regions overlap");
      XASSERTM((src < dst) || ((const char*)dst + bytes <= src), "source and destination memory regions overlap");

      ::memcpy(dst, src, bytes);
    }
  } // namespace Memory
} // namespace FEAT
