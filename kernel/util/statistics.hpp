#pragma once
#ifndef KERNEL_STATISTICS_HPP
#define KERNEL_STATISTICS_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>

namespace FEAST
{
  /**
   * \brief Statistics collection class
   *
   * The class Statistics encapsulates various hardware counters,
   * that collects e.g. the flop count of an linear solver.
   */
  class Statistics
  {
    private:

      /// global flop counter
      static Index _flops;

    public:

      /// Add an amount of flops to the global flop counter
      static void add_flops(Index flops);

      /// Retrieve current global flop counter
      static Index get_flops();

      /// Retrieve formated flops per second string
      static String get_formated_flops(double seconds);

      /// Reset global flop counter
      static void reset_flops();

  };
} // namespace FEAST

#endif // KERNEL_STATISTICS_HPP
