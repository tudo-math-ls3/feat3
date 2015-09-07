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
   * that collects e.g. the accumulated flop count of an linear solver.
   */
  class Statistics
  {
    private:

      /// global flop counter
      static Index _flops;

      /// global time of execution for reduction type operations
      static double _time_reduction;

      /// global time of execution for blas-2 type operations
      static double _time_spmv;

      /// global time of execution for blas-1 type operations
      static double _time_axpy;

      /// global time of execution for special preconditioner kernel type operations
      static double _time_precon;

      /// global time of execution for mpi related operations, e.g. send/recv or gather/scatter
      static double _time_mpi_execute;

      /// global time of execution for mpi related idle/wait tasks
      static double _time_mpi_wait;

    public:

      /// Add an amount of flops to the global flop counter
      static void add_flops(Index flops);

      /// Retrieve current global flop counter
      static Index get_flops();

      /// Retrieve formated flops per second string
      static String get_formated_flops(double seconds);

      /// Reset global flop counter
      static void reset_flops();

      static void add_time_reduction(double seconds);
      static void add_time_spmv(double seconds);
      static void add_time_axpy(double seconds);
      static void add_time_precon(double seconds);
      static void add_time_mpi_execute(double seconds);
      static void add_time_mpi_wait(double seconds);

      static double add_time_reduction();
      static double add_time_spmv();
      static double add_time_axpy();
      static double add_time_precon();
      static double add_time_mpi_execute();
      static double add_time_mpi_wait();

      /// Retrieve formated time consumption overview in percent
      static String get_formated_times();

      /// Retrieve formated time consumption overview in percent relative to some provided total time
      static String get_formated_times(double total_time);

      /// Reset all global timer counters
      static void reset_times();

  };
} // namespace FEAST

#endif // KERNEL_STATISTICS_HPP
