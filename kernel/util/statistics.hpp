#pragma once
#ifndef KERNEL_UTIL_STATISTICS_HPP
#define KERNEL_UTIL_STATISTICS_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/kahan_summation.hpp>

#include <list>
#include <map>
#include <time.h>
#include <iostream>
#include <fstream>

namespace FEAST
{
  /// Defect norm and time of execution for all iterations of a solver execution
  struct SolverStatistics
  {
    ///contains def_init and def_cur per iteration, implicitly contains the iteration count, too.
    /// \todo use proper datatype via template parameter, eg quad for quad based solvers
    std::vector<double> defect;

    /// Time of each iteration in seconds
    std::vector<double> toe;

    /// Time of mpi execution of each iteration in seconds
    std::vector<double> mpi_toe;
  };

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

      /// global time of execution for mpi related idle/wait tasks
      static KahanAccumulation _time_mpi_wait;

      /// map of SolverStatistics and their corresponding solver name
      static std::map<FEAST::String, SolverStatistics> _solver_statistics;

      static String _format_solver_statistics(String branch, SolverStatistics & stat)
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

        result +="MPI Iteration Timings [s]:\n";
        for (Index i(0) ; i < stat.mpi_toe.size() ; ++i)
        {
          if (i == 0 && stat.mpi_toe.at(i) == double(-1))
            continue;
          if (stat.mpi_toe.at(i) == double(-1))
            result += "------\n";
          else
            result += stringify_fp_sci(stat.mpi_toe.at(i)) + "\n";
        }

        return result;
      }

    public:

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

      /// Retrieve formated flops per second string
      static String get_formated_flops(double seconds, Index ranks = 1)
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
      inline static void add_time_mpi_wait(double seconds)
      {
        _time_mpi_wait = KahanSum(_time_mpi_wait, seconds);
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
      inline static double get_time_mpi_wait()
      {
        return _time_mpi_wait.sum;
      }

      /// add toe statistics entry for specific solver (branch name)
      inline static void add_solver_toe(String solver, double seconds)
      {
        auto it = _solver_statistics.find(solver);

        if (it != _solver_statistics.end())
        {
          it->second.toe.push_back(seconds);
        }
        else
        {
          SolverStatistics temp;
          temp.toe.push_back(seconds);
          _solver_statistics[solver] = temp;
        }
      }

      /// add mpi toe statistics entry for specific solver (branch name)
      inline static void add_solver_mpi_toe(String solver, double seconds)
      {
        auto it = _solver_statistics.find(solver);

        if (it != _solver_statistics.end())
        {
          it->second.mpi_toe.push_back(seconds);
        }
        else
        {
          SolverStatistics temp;
          temp.mpi_toe.push_back(seconds);
          _solver_statistics[solver] = temp;
        }
      }

      /// add defect norm statistics entry for specific solver (branch name)
      inline static void add_solver_defect(String solver, double defect)
      {
        auto it = _solver_statistics.find(solver);
        if (it != _solver_statistics.end())
        {
          it->second.defect.push_back(defect);
        }
        else
        {
          SolverStatistics temp;
          temp.defect.push_back(defect);
          _solver_statistics[solver] = temp;
        }
      }

      ///reset solver statistics (toe / defect norm)
      inline static void reset_solver_statistics()
      {
        _solver_statistics.clear();
      }

      /**
       * \brief Generate output with detailed solver statistics
       *
       * The generated String contains detailed detailed defect norms and execution times
       * for each used solver.
       *
       * The output per solver contains the following informations:
       * - solver name (branch name)
       * - defect norm of each iteration
       * - time of exection of each iteration
       *
       * \note All FEAST::Solver::IterativeSolver instances contain the initial defect norm as the first iteration defect norm.
       *
       * \note The FEAST::Solver::BasicVCycle prints its TOE from coarse to fine level.
       */
      static String get_formated_solvers()
      {
        String result;

        for (auto stat : _solver_statistics)
        {
          result += _format_solver_statistics(stat.first, stat.second);
          result += "\n";
        }

        return result;
      }

      /**
       * \copydoc get_formated_solvers()
       *
       * \param[in] branch The solver (as complete branch name) that shall be used as the tree's root.
       */
      static String get_formated_solver(String branch)
      {
        auto it = _solver_statistics.find(branch);
        if (it == _solver_statistics.end())
          throw InternalError("solver with branch " + branch + " not found!");
        else
          return _format_solver_statistics(branch, it->second);
      }

      /*
      /// Retrieve formated time consumption overview in percent
      static String get_formated_times()
      {
        double total_time = get_time_reduction() + get_time_spmv() + get_time_axpy() + get_time_precon() + get_time_mpi_execute() + get_time_mpi_wait();
        return get_formated_times(total_time);
      }*/

      /// Retrieve formated time consumption overview in percent relative to some provided total time
      static String get_formated_times(double total_time)
      {
        String result = "Total time: " + stringify(total_time) + "s";
        if (total_time == 0.)
          return result;

        KahanAccumulation measured_time;
        measured_time = KahanSum(measured_time, get_time_reduction());
        measured_time = KahanSum(measured_time, get_time_spmv());
        measured_time = KahanSum(measured_time, get_time_axpy());
        measured_time = KahanSum(measured_time, get_time_precon());
        measured_time = KahanSum(measured_time, get_time_mpi_execute());
        /// \todo readd mpi_wait time, once wait is no more included in execute
        // mpi_wait is included in mpi_execute
        //measured_time = KahanSum(measured_time, get_time_mpi_wait());
        // subtract wait from execute time, as the later is included in the former accidentaly
        _time_mpi_execute.sum -= _time_mpi_wait.sum;

        if (measured_time.sum > total_time)
          throw InternalError("Accumulated op time (" + stringify(measured_time.sum) + ") is greater as the provided total execution time (" + stringify(total_time) + ") !");

        result += "\n";
        result += "Accumulated op time: " + stringify(measured_time.sum) + "\n";

        result += "\n";
        result += String("Reductions:").pad_back(17) + stringify(get_time_reduction() / total_time * 100.) + "%\n";
        result += String("Blas-1:").pad_back(17) + stringify(get_time_axpy() / total_time * 100.) + "%\n";
        result += String("Blas-2:").pad_back(17) + stringify(get_time_spmv() / total_time * 100.) + "%\n";
        result += String("Precon Kernels:").pad_back(17) + stringify(get_time_precon() / total_time * 100.) + "%\n";
        result += String("MPI Execution:").pad_back(17) + stringify(get_time_mpi_execute() / total_time * 100.) + "%\n";
        result += String("MPI Wait:").pad_back(17) + stringify(get_time_mpi_wait() / total_time * 100.) + "%\n";
        result += String("Not covered:").pad_back(17) + stringify( (total_time - measured_time.sum) / total_time * 100.) + "%";
        return result;
      }

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
        _time_mpi_wait.sum = 0.;
        _time_mpi_wait.correction = 0.;
      }

      static void write_out_solver_statistics(Index rank, size_t la_bytes, size_t domain_bytes, size_t mpi_bytes, String filename = "solver_stats")
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
        file << "reduction " << stringify(get_time_reduction()) << std::endl;
        file << "axpy " << stringify(get_time_axpy()) << std::endl;
        file << "spmv " << stringify(get_time_spmv()) << std::endl;
        file << "precon " << stringify(get_time_precon()) << std::endl;
        file << "mpi " << stringify(get_time_mpi_execute()) << std::endl;

        file << std::endl;

        file << "#global size" << std::endl;
        file << "la " << stringify(la_bytes) << std::endl;
        file << "domain " << stringify(domain_bytes) << std::endl;
        file << "mpi " << stringify(mpi_bytes) << std::endl;

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

          file << "#mpi_toe" << std::endl;
          for (Index i(0) ; i < stat.second.mpi_toe.size() ; ++i)
          {
            if (i == 0 && stat.second.mpi_toe.at(i) == double(-1))
              continue;
            if (stat.second.mpi_toe.at(i) == double(-1))
              file << "-" << std::endl;
            else
              file << stringify_fp_sci(stat.second.mpi_toe.at(i)) << std::endl;
          }
        }

        file.close();
      }
  };
} // namespace FEAST

#endif // KERNEL_UTIL_STATISTICS_HPP
