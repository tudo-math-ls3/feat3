#pragma once
#ifndef KERNEL_LAFEM_ABSOLUTE_HPP
#define KERNEL_LAFEM_ABSOLUTE_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>

#include <cmath>


namespace FEAST
{
  namespace LAFEM
  {
    template<typename DT_>
    struct Absolute
    {
      static DT_ value(DT_ val)
      {
        return std::abs(val);
      }
    };

    template<>
    struct Absolute<Index>
    {
      static Index value(Index val)
      {
        return val;
      }

    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ABSOLUTE_HPP
