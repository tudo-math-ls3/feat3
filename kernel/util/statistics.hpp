#pragma once
#ifndef KERNEL_UTIL_STATISTICS_HPP
#define KERNEL_UTIL_STATISTICS_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/kahan_summation.hpp>
#include <kernel/util/comm_base.hpp>
#include <kernel/solver/expression.hpp>

#include <list>
#include <map>
#include <time.h>
#include <iostream>
#include <fstream>
#include <memory>

namespace FEAT
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
      static KahanAccumulation _time_reduction;

      /// global time of execution for blas-2 type operations
      static KahanAccumulation _time_spmv;

      /// global time of execution for blas-1 type operations
      static KahanAccumulation _time_axpy;

      /// global time of execution for special preconditioner kernel type operations
      static KahanAccumulation _time_precon;

      /// global time of execution for mpi related operations, e.g. send/recv or gather/scatter
      static KahanAccumulation _time_mpi_execute;

      /// global time of execution for mpi related idle/wait tasks of reduction operations
      static KahanAccumulation _time_mpi_wait_reduction;

      /// global time of execution for mpi related idle/wait tasks of spmv operations
      static KahanAccumulation _time_mpi_wait_spmv;

      /// a consecutive list of all solver actions
      static std::map<String, std::list<std::shared_ptr<Solver::ExpressionBase>>> _solver_expressions;

      /*static String _format_solver_statistics(String branch, SolverStatistics & stat)
      {
        String result;

        result += "\n" + branch + "\n";
        result +="Iteration Defects:\n";

        for (Index i(0) ; i < stat.defect.size() ; ++i)
        {
          if (i == 0 && stat.defect.at(i) == double(-1))
            continue;
          if (stat.defect.at(i) == double(-1))
            result += "------\n";
          else
            result += stringify_fp_sci(stat.defect.at(i)) + "\n";
        }

        result +="Iteration Timings [s]:\n";
        for (Index i(0) ; i < stat.toe.size() ; ++i)
        {
          if (i == 0 && stat.toe.at(i) == double(-1))
            continue;
          if (stat.toe.at(i) == double(-1))
            result += "------\n";
          else
            result += stringify_fp_sci(stat.toe.at(i)) + "\n";
        }

        result +="MPI Execution Iteration Timings [s]:\n";
        for (Index i(0) ; i < stat.mpi_execute.size() ; ++i)
        {
          if (i == 0 && stat.mpi_execute.at(i) == double(-1))
            continue;
          if (stat.mpi_execute.at(i) == double(-1))
            result += "------\n";
          else
            result += stringify_fp_sci(stat.mpi_execute.at(i)) + "\n";
        }

        result +="MPI Wait Iteration Timings [s]:\n";
        for (Index i(0) ; i < stat.mpi_wait.size() ; ++i)
        {
          if (i == 0 && stat.mpi_wait.at(i) == double(-1))
            continue;
          if (stat.mpi_wait.at(i) == double(-1))
            result += "------\n";
          else
            result += stringify_fp_sci(stat.mpi_wait.at(i)) + "\n";
        }

        return result;
      }*/

    public:

      /// time of partitioning in seconds, needs initialisation
      static double toe_partition;
      /// time of assembly in seconds, needs initialisation
      static double toe_assembly;
      /// time of solution in seconds, needs initialisation
      static double toe_solve;
      /// the current solver's descriptive string
      static String expression_target;


      /// Reset all internal counters (flops/times/solver_statistics)
      static void reset()
      {
        reset_flops();
        reset_times();
        reset_solver_statistics();
      }

      /// Add an amount of flops to the global flop counter
      static void add_flops(Index flops)
      {
        _flops += flops;
      }

      /// Retrieve current global flop counter
      static Index get_flops()
      {
        return _flops;
      }

      /// Retrieve formatted flops per second string
      static String get_formatted_flops(double seconds, Index ranks = 1)
      {
        double flops((double)_flops);
        flops /= seconds;
        flops /= 1000.; // kilo
        flops /= 1000.; // mega
        flops /= 1000.; // giga
        return stringify(flops) + " GFlop/s [" + stringify(flops * double(ranks)) + " GFlops/s]";
      }

      /// Reset global flop counter
      static void reset_flops()
      {
        _flops = Index(0);
      }

      inline static void add_time_reduction(double seconds)
      {
        _time_reduction = KahanSum(_time_reduction, seconds);
      }
      inline static void add_time_spmv(double seconds)
      {
        _time_spmv = KahanSum(_time_spmv, seconds);
      }
      inline static void add_time_axpy(double seconds)
      {
        _time_axpy = KahanSum(_time_axpy, seconds);
      }
      inline static void add_time_precon(double seconds)
      {
        _time_precon = KahanSum(_time_precon, seconds);
      }
      inline static void add_time_mpi_execute(double seconds)
      {
        _time_mpi_execute = KahanSum(_time_mpi_execute, seconds);
      }
      inline static void add_time_mpi_wait_reduction(double seconds)
      {
        _time_mpi_wait_reduction = KahanSum(_time_mpi_wait_reduction, seconds);
      }
      inline static void add_time_mpi_wait_spmv(double seconds)
      {
        _time_mpi_wait_spmv = KahanSum(_time_mpi_wait_spmv, seconds);
      }

      inline static double get_time_reduction()
      {
        return _time_reduction.sum;
      }
      inline static double get_time_spmv()
      {
        return _time_spmv.sum;
      }
      inline static double get_time_axpy()
      {
        return _time_axpy.sum;
      }
      inline static double get_time_precon()
      {
        return _time_precon.sum;
      }
      inline static double get_time_mpi_execute()
      {
        return _time_mpi_execute.sum;
      }
      inline static double get_time_mpi_wait_reduction()
      {
        return _time_mpi_wait_reduction.sum;
      }
      inline static double get_time_mpi_wait_spmv()
      {
        return _time_mpi_wait_spmv.sum;
      }

      inline static void add_solver_expression(std::shared_ptr<Solver::ExpressionBase> expression)
      {
        _solver_expressions[expression_target].push_back(expression);
      }

      /**
       * \brief Returns a descriptive string of the complete solver tree.
       *
       * Create and format a string describing the complete solver tree.
       *
       * \note This method makes some simplifications, e.g. stating only one smoother
       * for the complete FEAT::Solver::BasicVCycle.
       *
       * \note The solver must have been executed (successfully) at least one time after the last reset_solver_statistics() call.
       */
      static String get_formatted_solver_tree();

      /// print out the complete solver expression list
      static void print_solver_expressions();

      ///reset solver statistics (toe / defect norm)
      inline static void reset_solver_statistics()
      {
        _solver_expressions.clear();
      }

      /*
      /// Retrieve formatted time consumption overview in percent
      static String get_formatted_times()
      {
        double total_time = get_time_reduction() + get_time_spmv() + get_time_axpy() + get_time_precon() + get_time_mpi_execute() + get_time_mpi_wait();
        return get_formatted_times(total_time);
      }*/

      /**
       * \brief Retrieve formatted time consumption overview in percent relative to some provided total time
       *
       * \note This method uses mpi collectives and thus needs to be called by all ranks, even if you dont use the result on every rank on your own.
       */
      static String get_formatted_times(double total_time);

      /// Reset all global timer counters
      static void reset_times()
      {
        _time_reduction.sum = 0.;
        _time_reduction.correction = 0.;
        _time_spmv.sum = 0.;
        _time_spmv.correction = 0.;
        _time_axpy.sum = 0.;
        _time_axpy.correction = 0.;
        _time_precon.sum = 0.;
        _time_precon.correction = 0.;
        _time_mpi_execute.sum = 0.;
        _time_mpi_execute.correction = 0.;
        _time_mpi_wait_reduction.sum = 0.;
        _time_mpi_wait_reduction.correction = 0.;
        _time_mpi_wait_spmv.sum = 0.;
        _time_mpi_wait_spmv.correction = 0.;
      }

      /*
      static void write_out_solver_statistics_scheduled(Index rank, size_t la_bytes, size_t domain_bytes, size_t mpi_bytes, Index cells, Index dofs, Index nzes, String filename = "solver_stats");

      static void write_out_solver_statistics(Index rank, size_t la_bytes, size_t domain_bytes, size_t mpi_bytes, Index cells, Index dofs, Index nzes, String filename = "solver_stats")
      {
        filename += ".";
        filename += stringify(rank);

        std::ofstream file(filename.c_str(), std::ofstream::out);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open statistics file " + filename);

        {
          time_t t = time(nullptr);
          struct tm * ts = localtime(&t);

          file << "timestamp " << stringify(ts->tm_mon + 1) << "/" << stringify(ts->tm_mday) << "/" << stringify(ts->tm_year + 1900) << " " << stringify(ts->tm_hour)
            << ":" << stringify(ts->tm_min) << ":" << stringify(ts->tm_sec) << std::endl;
        }

        file << std::endl;

        file << "#global time" << std::endl;
        file << "partition " << stringify(toe_partition) << std::endl;
        file << "assembly " << stringify(toe_assembly) << std::endl;
        file << "solve " << stringify(toe_solve) << std::endl;

        file << std::endl;

        file << "#solver time" << std::endl;
        file << "reduction " << stringify(get_time_reduction()) << std::endl;
        file << "axpy " << stringify(get_time_axpy()) << std::endl;
        file << "spmv " << stringify(get_time_spmv()) << std::endl;
        file << "precon " << stringify(get_time_precon()) << std::endl;
        file << "mpi " << stringify(get_time_mpi_execute()) << std::endl;
        file << "wait reduction" << stringify(get_time_mpi_wait_reduction()) << std::endl;
        file << "wait spmv" << stringify(get_time_mpi_wait_spmv()) << std::endl;

        file << std::endl;

        file << "#global size" << std::endl;
        file << "domain " << stringify(domain_bytes) << std::endl;
        file << "mpi " << stringify(mpi_bytes) << std::endl;
        file << "la " << stringify(la_bytes) << std::endl;

        file << std::endl;

        file <<"#problem size" << std::endl;
        file << "cells " << stringify(cells) << std::endl;
        file << "dofs " << stringify(dofs) << std::endl;
        file << "nzes " << stringify(nzes) << std::endl;

        file << std::endl;

        file << "#solver statistics" << std::endl;
        for (auto stat : _solver_statistics)
        {
          file << stat.first << std::endl;

          file << "#defects" << std::endl;
          for (Index i(0) ; i < stat.second.defect.size() ; ++i)
          {
            if (i == 0 && stat.second.defect.at(i) == double(-1))
              continue;
            if (stat.second.defect.at(i) == double(-1))
              file << "-" << std::endl;
            else
              file << stringify_fp_sci(stat.second.defect.at(i)) << std::endl;
          }

          file << "#toe" << std::endl;
          for (Index i(0) ; i < stat.second.toe.size() ; ++i)
          {
            if (i == 0 && stat.second.toe.at(i) == double(-1))
              continue;
            if (stat.second.toe.at(i) == double(-1))
              file << "-" << std::endl;
            else
              file << stringify_fp_sci(stat.second.toe.at(i)) << std::endl;
          }

          file << "#mpi_execute" << std::endl;
          for (Index i(0) ; i < stat.second.mpi_execute.size() ; ++i)
          {
            if (i == 0 && stat.second.mpi_execute.at(i) == double(-1))
              continue;
            if (stat.second.mpi_execute.at(i) == double(-1))
              file << "-" << std::endl;
            else
              file << stringify_fp_sci(stat.second.mpi_execute.at(i)) << std::endl;
          }

          file << "#mpi_wait" << std::endl;
          for (Index i(0) ; i < stat.second.mpi_wait.size() ; ++i)
          {
            if (i == 0 && stat.second.mpi_wait.at(i) == double(-1))
              continue;
            if (stat.second.mpi_wait.at(i) == double(-1))
              file << "-" << std::endl;
            else
              file << stringify_fp_sci(stat.second.mpi_wait.at(i)) << std::endl;
          }
        }

        file.close();
      }*/
  };
} // namespace FEAT

#endif // KERNEL_UTIL_STATISTICS_HPP
