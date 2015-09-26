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

      inline static void add_time_reduction(double seconds)
      {
        _time_reduction += seconds;
      }
      inline static void add_time_spmv(double seconds)
      {
        _time_spmv += seconds;
      }
      inline static void add_time_axpy(double seconds)
      {
        _time_axpy += seconds;
      }
      inline static void add_time_precon(double seconds)
      {
        _time_precon += seconds;
      }
      inline static void add_time_mpi_execute(double seconds)
      {
        _time_mpi_execute += seconds;
      }
      inline static void add_time_mpi_wait(double seconds)
      {
        _time_mpi_wait += seconds;
      }

      inline static double get_time_reduction()
      {
        return _time_reduction;
      }
      inline static double get_time_spmv()
      {
        return _time_spmv;
      }
      inline static double get_time_axpy()
      {
        return _time_axpy;
      }
      inline static double get_time_precon()
      {
        return _time_precon;
      }
      inline static double get_time_mpi_execute()
      {
        return _time_mpi_execute;
      }
      inline static double get_time_mpi_wait()
      {
        return _time_mpi_wait;
      }

      /// Retrieve formated time consumption overview in percent
      static String get_formated_times();

      /// Retrieve formated time consumption overview in percent relative to some provided total time
      static String get_formated_times(double total_time);

      /// Reset all global timer counters
      static void reset_times();

  };
} // namespace FEAST

#endif // KERNEL_STATISTICS_HPP
